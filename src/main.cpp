#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

const double velocity_target = 10.;

// for convenience
using json = nlohmann::json;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

double prev_x, prev_y, prev_v;
chrono::system_clock::time_point prev_t = chrono::system_clock::now();

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;
  Eigen::VectorXd coeffs(3);

//  coeffs << 5, 0, 0;
//  mpc.Solve(coeffs, 0, 10);

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
//    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {

          // Parse data

          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = (double)j[1]["speed"] * 0.44704; // Convert speed to m/s

          // Measure velocity and acceleration (for debugging and inspection)

          chrono::system_clock::time_point t = chrono::system_clock::now();
          double dtime = chrono::duration_cast<chrono::milliseconds>(t - prev_t).count() / 1000.0;
          auto ddistance = sqrt(pow(px - prev_x, 2) + pow(py - prev_y, 2));
          auto dvelocity = v - prev_v;
          cout << "v=" << v << ", ds/dt=" << ddistance/dtime << ", dv/dt=" << dvelocity/dtime << ", dt=" << dtime << endl;
          prev_x = px;
          prev_y = py;
          prev_v = v;
          prev_t = t;

          // Convert to car-coordinate system

          vector<double> relx;
          vector<double> rely;

          for (size_t i = 0; i < ptsx.size(); i++) {
            double x = ptsx.at(i) - px;
            double y = ptsy.at(i) - py;
            double xrot = x * cos(psi) + y * sin(psi);
            double yrot = y * cos(psi) - x * sin(psi);
            relx.push_back(xrot);
            rely.push_back(yrot);
          }

          // Polynomial fit

          Eigen::Map<Eigen::VectorXd> relx_vec(relx.data(), relx.size());
          Eigen::Map<Eigen::VectorXd> rely_vec(rely.data(), rely.size());
          auto coeffs = polyfit(relx_vec, rely_vec, 3);

          // MPC solve

          double delta, a;
          vector<double> mpc_x, mpc_y;
          tie(delta, a, mpc_x, mpc_y) = mpc.Solve(coeffs, v, velocity_target);

          // Map (delta,a) to (steer,throttle)

          // Steering angle is negative of delta (1. == -25 degrees)
          double steer_value = - delta / (25 / 180.0 * M_PI);
          // TODO: not quite correct, throttle is not proportional to acceleration
          // throttle=1 at velocity=0 results in 5m/s2 acceleration
          double throttle_value = a / 5;

          // Steer and Throttle must be in [-1, 1]
          if (steer_value < -1) steer_value = -1;
          if (steer_value > 1) steer_value = 1;
          if (throttle_value < -1) throttle_value = -1;
          if (throttle_value > 1) throttle_value = 1;
          cout << "steering=" << steer_value << ", throttle=" << throttle_value << endl;

          // Form message

          json msgJson;
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;
          msgJson["mpc_x"] = mpc_x;
          msgJson["mpc_y"] = mpc_y;
          msgJson["next_x"] = relx;
          msgJson["next_y"] = rely;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
//          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
