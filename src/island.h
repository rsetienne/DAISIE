#pragma once

struct island {
  virtual void update_area(double t) = 0;
  double area;
};

struct island_static : island {
  void update_area(double t) override {
    // static nothing happens
    return;
  }
  
  double area = 1.0;
};

class A_beta {
public:
  A_beta(double proptime_max,
         double max_time,
         double max_area) :
    f(proptime_max / (1 - proptime_max)),
      Tmax(max_time),
      Amax(max_area) {
  }
  
  double calc_Abeta(double t, double peak) const {
    double proptime = t / Tmax;
    double a =  f * peak / (1 + f);
    double b = peak / (1 + f); 
    
    double divisor = std::pow(a / (a + b), a) * std::pow((b / (a + b)), b);
    double x = std::pow(1 - proptime, b) / divisor;
    return Amax * std::pow(proptime, a) * x;
  }
  
private:
  const double f;
  const double Tmax;
  const double Amax;
};

class angular {
public:
  angular(double sea_level_frequency,
          double total_island_age,
          double total_time,
          double amplitude,
          double current_area,
          double island_gradient_angle) : 
    angular_freq(2 * M_PI * sea_level_frequency),
    proptime_curr(total_time / total_island_age),
    ampl(amplitude),
    theta(island_gradient_angle),
    current_area_(current_area) {
  }

  double calc_angular(double proptime, double Acurr = -1) const {
    double Acurr_ = Acurr < 0.0 ? current_area_ : Acurr;
    auto delta_sl = ampl * std::cos((proptime_curr - proptime) * angular_freq);
    auto r_curr = std::sqrt((Acurr_) / M_PI);
    auto h_curr = std::tan(theta) * r_curr;
    auto h_delta = h_curr - ampl + delta_sl;
    if (h_delta < 0) h_delta = 0.0;
    auto x = h_delta / std::tan(theta);
    return M_PI * x * x;
  }
  
private:
  const double angular_freq;
  const double proptime_curr;
  const double ampl;
  const double theta;
  const double current_area_;
};


struct island_beta : island{
public:
  island_beta(double max_area,
              double total_island_age,
              double proportional_peak_t) : 
              Abeta_(proportional_peak_t, total_island_age, max_area) {
  }
  
  double area = 1.0;
  
  void update_area(double t) override {
    // TODO: make peak correct
    double peak = 1.0;
    area = Abeta_.calc_Abeta(t, peak);
  }
private:
  A_beta Abeta_;
};

struct island_angular : island {
public:
  island_angular(double sea_level_frequency,
                 double total_island_age,
                 double total_time,
                 double amplitude,
                 double current_area,
                 double island_gradient_angle) : 
  Tmax(total_island_age),
  Acurr(current_area),
    angular_(angular(sea_level_frequency,
                     total_island_age,
                     total_time,
                     amplitude,
                     current_area,
                     island_gradient_angle))
    {
  }
  
  double area = 1.0;
  
  void update_area(double t) override {
    double proptime = t / Tmax;
    area = angular_.calc_angular(proptime, Acurr);
  }
private:
  const double Tmax;
  const double Acurr;
  angular angular_;
};


struct island_beta_angular : island {
public:
  island_beta_angular(double sea_level_frequency,
                      double total_island_age,
                      double total_time,
                      double amplitude,
                      double current_area,
                      double island_gradient_angle,
                      double max_area,
                      double proportional_peak_t) : 
  angular_(angular(sea_level_frequency,
                   total_island_age,
                   total_time,
                   amplitude,
                   current_area,
                   island_gradient_angle)),
    abeta_(A_beta(proportional_peak_t, total_island_age, max_area)),
    Tmax(total_time) {
  }
  
  double area = 1.0;
  
  void update_area(double t) override {
    double peak = 1.0;
    double proptime = t / Tmax;
    auto abeta_val = abeta_.calc_Abeta(t, peak);
    area = angular_.calc_angular(proptime, abeta_val);
  }
  
private:
  angular angular_;
  A_beta abeta_;
  const double Tmax;
};

