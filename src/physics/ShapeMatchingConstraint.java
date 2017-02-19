package physics;

import java.util.LinkedList;

public class ShapeMatchingConstraint {
    public int id;
    public LinkedList<Particle> particles = new LinkedList<>();
    public double stiffness;
    public double centerOfMass[] = {0, 0};

    public ShapeMatchingConstraint(int id, Particle pa[], double stiffness) {
        this.id = id;
        this.stiffness = stiffness;
        for (Particle p : pa) {
            particles.add(p);
            centerOfMass[0] += p.x;
            centerOfMass[1] += p.y;
        }
        centerOfMass[0] /= particles.size();
        centerOfMass[1] /= particles.size();
    }
}
