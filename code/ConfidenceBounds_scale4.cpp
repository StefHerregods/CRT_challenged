// To use the simulation function, install the Rcpp library in R and source this file with - sourceCpp("code/ConfidenceBounds_scale4.cpp") -

// Defining libraries, dependencies and make function available for use in R
#include <Rcpp.h>
// [[Rcpp::depends(RcppZiggurat)]]
#include <Ziggurat.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;
// [[Rcpp::export]]

// Simulation function
NumericMatrix DDM_confidence_bounds(int ntrials, double s, double dt, double z, double a, double v, double a_slope, double ter, double a2, double postdriftmod, double a2_slope_upper, double a2_slope_lower, double ter2, double starting_point_confidence) {
  
  // Initialize output
  NumericMatrix DATA(ntrials,4);
  
  // Loop over trials
  for (int i = 0; i < ntrials; i++) {

    // Initialize variables
    double evidence = a*z;
    double t = 0;
    double t2 = 0;
    int cor = -1;

    // Decisional processing
    while (true){
      
      // Accumulate evidence over time
      t = t + dt;
      evidence = evidence + v * dt + s * sqrt(dt) * zigg.norm();
      
      // Make a decision once a boundary is crossed
      if (evidence >= a - t * a_slope){
        cor = 1;
        break;
      } else if (evidence <= 0 + t * a_slope) {
        cor = 0;
        break;
      }
    }

    // Save decision data
    DATA(i,0) = (t + ter);
    DATA(i,1) = cor;

    // Restart evidence accumulation
    double v_post = v * postdriftmod;
    evidence = a2 * starting_point_confidence;
    if (cor == 0){
      v_post = -1 * v_post;
    }

    // Post-decisional evidence accumulation until reaching a confidence boundary
    while ((evidence < a2 - t2*a2_slope_upper) && (evidence > t2*a2_slope_lower)){
      t2 = t2 + dt;
      evidence = evidence + v_post * dt + s * sqrt(dt) * zigg.norm();
    }

    // Save confidence RT
    DATA(i,2) = t2 + ter2;

    // Save confidence judgments
    if (evidence < a2/4){
      DATA(i,3) = 1;
    } else if (evidence < 2*a2/4){
      DATA(i,3) = 2;
    } else if (evidence < 3*a2/4){
      DATA(i,3) = 3;
    } else {
      DATA(i,3) = 4;
    }

  }
  return DATA; // decision RT, decision accuracy, confidence RT, confidence judgment
}