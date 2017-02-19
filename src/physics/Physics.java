package physics;

import java.awt.BorderLayout;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.AffineTransform;
import javax.swing.JFrame;
import javax.swing.JPanel;

// http://labs.byhook.com/2010/06/29/particle-based-rigid-bodies-using-shape-matching/
// https://github.com/ClickerMonkey/ImpulseEngine/blob/master/src/org/magnos/impulse/CollisionPolygonPolygon.java
// http://research.nii.ac.jp/~takayama/teaching/utokyo-iscg-2015/assignments/sample/iscg-2015-assignment-c-pbd2d-sample.html
public class Physics extends JPanel {

    private static final double dt = .05d;
    private static final Scene scene = new Scene();

    public static void main(String[] args) {

        JFrame frame = new JFrame("HoE Physics");
        frame.setSize(1200, 800);
        frame.setLocationRelativeTo(null);
        frame.setLayout(new BorderLayout());
        frame.add(new Physics(), BorderLayout.CENTER);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
        
        // Jelenet létrehozása.
        double res = .5;
        double fric = .1;
        scene.addParticle(2, 10, 300, res, fric);
        scene.addParticle(1, 80, 300, res, fric);
        scene.addWall(0, 0, 0, 1);
        scene.addWall(0, 0, Math.cos(2.8), Math.sin(2.8));
        scene.addWall(0, 50, Math.cos(2), Math.sin(2));

        // Szimuláció futtatása.
        new Thread(() -> {
            while (true) {
                // Léptetés.
                scene.step(dt);
                
                // Kirajzolása.
                frame.repaint();
                try {
                    Thread.sleep((long) (dt * 100));
                } catch (InterruptedException ex) {}
            }
        }).start();

    }

    @Override
    protected void paintComponent(Graphics g) {

        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g.create();
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        AffineTransform currentTransform = g2.getTransform();
        AffineTransform viewTransform = new AffineTransform(currentTransform);
        viewTransform.scale(1, -1);
        viewTransform.translate(getWidth() / 2, -getHeight() * 7 / 8);
        g2.setTransform(viewTransform);

        // Jelenet kirajzolása.
        scene.draw(g2);
    }

}
