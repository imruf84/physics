package physics;

public class Particle {
    public int id;
    public double mass;
    public double x;
    public double y;
    public double tempX;
    public double tempY;
    public double vx;
    public double vy;
    public double restitution = .05;
    public double friction = .5;

    public Particle(int id, double mass, double x, double y, double vx, double vy, double restitution, double friction) {
        this.id = id;
        this.mass = mass;
        this.x = x;
        this.y = y;
        this.vx = vx;
        this.vy = vy;
        this.restitution = restitution;
        this.friction = friction;
    }

    
}
