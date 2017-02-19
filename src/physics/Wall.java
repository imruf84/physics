package physics;

public class Wall {
    public int id;
    public double x;
    public double y;
    public double nx;
    public double ny;

    public Wall(int id, double x, double y, double nx, double ny) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.nx = nx;
        this.ny = ny;
    }
    
    public double dist(double px, double py) {
        return nx * (px - x) + ny * (py - y);
    }
}
