class Particle {
  private:
    double mass;
    double position[3], velocity[3];
  public:
    Particle(double imass=0.0) {
      mass = imass;
      for (int i=0; i<3; i++) {
        position[i] = 0.0;
        velocity[i] = 0.0;
      }
    }
    virtual ~Particle() { }
    virtual double Position(int i) const { return position[i]; }
    virtual double Velocity(int i) const { return velocity[i]; }
    virtual double Kinetic_Energy() const {
      double ke = 0.0;
      for (int i=0; i<3; i++) ke += velocity[i]*velocity[i];
      return mass*ke;
    }
    virtual double Charge() const { return 0.0; }
};
