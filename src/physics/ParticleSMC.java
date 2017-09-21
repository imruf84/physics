package physics;

public class ParticleSMC extends Particle {
    
    public double rxRest;
    public double ryRest;
    public double rx;
    public double ry;
    
    public ParticleSMC(int id, double mass, double x, double y, double vx, double vy, double restitution, double friction) {
        super(id, mass, x, y, vx, vy, restitution, friction);
    }
    
}
