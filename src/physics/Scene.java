package physics;

import java.awt.Graphics2D;
import java.util.LinkedList;

public class Scene {
    
    double gravity[] = {.0, -9.8};
    double damping = 0;
    
    int wallIdCounter = 0;
    LinkedList<ParticleSMC> particles = new LinkedList<>();
    LinkedList<Wall> walls = new LinkedList<>();
    LinkedList<Constraints> constraints = new LinkedList<>();
    
    Particle addParticle(double mass, double x, double y, double vx, double vy, double restitution, double friction) {
        
        ParticleSMC p = new ParticleSMC(particles.size(), mass, x, y, vx, vy, restitution, friction);
        particles.add(p);
        
        return p;
    }
    
    Wall addWall(double x, double y, double nx, double ny) {
        Wall w = new Wall(wallIdCounter++, x, y, nx, ny);
        walls.add(w);
        
        return w;
    }
    
    Constraints addConstraints(Constraints c) {
        constraints.add(c);
        
        return c;
    }
    
    public void step(double dt) {
        
        // Erők.
        for (Particle p : particles) {
            double force[] = {0, 0};
            force[0] = gravity[0];
            force[1] = gravity[1];
            
            p.vx += force[0] * dt / p.mass;
            p.vy += force[1] * dt / p.mass;
            
            p.vx *= 1 - damping;
            p.vy *= 1 - damping;
            
            p.tempX = p.x + dt * p.vx;
            p.tempY = p.y + dt * p.vy;
            
            // Ütközés a falakkal.
            for (Wall w : walls) {
                double dist = w.dist(p.tempX, p.tempY);
                if (dist <= 0) {
                    double vnx = w.nx * (w.nx * p.vx + w.ny * p.vy);
                    double vny = w.ny * (w.nx * p.vx + w.ny * p.vy);
                    double vtx = p.vx - vnx;
                    double vty = p.vy - vny;
                    vnx *= -p.restitution;
                    vny *= -p.restitution;
                    vtx *= 1 - p.friction;
                    vty *= 1 - p.friction;
                    p.tempX = p.x + dt * (vnx + vtx);
                    p.tempY = p.y + dt * (vny + vty);
                    dist = w.dist(p.tempX, p.tempY);
                    if (dist < 0) {
                        p.tempX += -dist * w.nx;
                        p.tempY += -dist * w.ny;
                    }
                }
            }
        }
        
        // Kényszerítők.
        for (int i = 0; i < Physics.solver_iterations; ++i) {
            for (Constraints c : constraints) {
                c.enforce();
            }
        }
        
        // Sebesség és pozíció frissítése.
        for (Particle p : particles) {
            p.vx = (p.tempX - p.x) * 1 / dt;
            p.vy = (p.tempY - p.y) * 1 / dt;
            p.x = p.tempX;
            p.y = p.tempY;
        }
    }
    
    public void draw(Graphics2D g2) {
        double particleRadius = 10;
        for (Particle p : particles) {
            g2.drawOval((int) (p.x - particleRadius / 2d), (int) (p.y - particleRadius / 2d), (int) particleRadius, (int) particleRadius);
        }
        
        for (Wall w : walls) {
            double d[] = {1000 * w.ny, -1000 * w.nx};
            g2.drawLine((int)(w.x+d[0]), (int)(w.y+d[1]), (int)(w.x-d[0]), (int)(w.y-d[1]));
        }
        
        for (Constraints c : constraints) {
            // TODO: szakaszok megjelenítése (síkidom oldalai)
        }
    }
    
    public ParticleSMC getParticle(int index) {
        return particles.get(index);
    }
}
