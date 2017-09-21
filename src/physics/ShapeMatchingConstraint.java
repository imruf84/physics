package physics;

import java.util.Arrays;
import java.util.LinkedList;

public class ShapeMatchingConstraint implements Constraints {

    public int id;
    public LinkedList<ParticleSMC> particles = new LinkedList<>();
    public double stiffness;
    public double[] centerOfMassRest = {0, 0};

    public ShapeMatchingConstraint(int id, ParticleSMC pa[], double stiffness) {
        this.id = id;
        this.stiffness = stiffness;
        for (ParticleSMC p : pa) {
            particles.add(p);
        }
        double com[] = calculateCenterOfMass();
        centerOfMassRest[0] = com[0];
        centerOfMassRest[1] = com[1];

        for (ParticleSMC p : particles) {
            p.rxRest = p.x - centerOfMassRest[0];
            p.ryRest = p.y - centerOfMassRest[1];
        }
    }

    public final double[] calculateCenterOfMass() {
        double com[] = {0, 0};
        for (ParticleSMC p : particles) {
            com[0] += p.x;
            com[1] += p.y;
        }

        com[0] /= particles.size();
        com[1] /= particles.size();

        return com;
    }

    public final double[] calculateCenterOfMassTemp() {
        double com[] = {0, 0};
        for (ParticleSMC p : particles) {
            com[0] += p.tempX;
            com[1] += p.tempY;
        }

        com[0] /= particles.size();
        com[1] /= particles.size();

        return com;
    }

    @Override
    public void enforce() {

        double com[] = calculateCenterOfMassTemp();
        double A[][] = {{0, 0}, {0, 0}};

        for (ParticleSMC p : particles) {
            p.rx = p.tempX - com[0];
            p.ry = p.tempY - com[1];

            A[0][0] += p.mass * p.rx * p.rxRest;
            A[0][1] += p.mass * p.rx * p.ryRest;
            A[1][0] += p.mass * p.ry * p.rxRest;
            A[1][1] += p.mass * p.ry * p.ryRest;
        }

        double AtA[][] = {
            {A[0][0] * A[0][0] + A[1][0] * A[1][0], A[0][0] * A[0][1] + A[1][0] * A[1][1]},
            {A[0][0] * A[0][1] + A[1][0] * A[1][1], A[0][1] * A[0][1] + A[1][1] * A[1][1]}};

        // SVD
        double temp;
//Compute the thin SVD from G. H. Golub and C. Reinsch, Numer. Math. 14, 403-420 (1970)
        double prec = Math.pow(2, -52); // assumes double prec
        double tolerance = 1.e-64 / prec;
        int itmax = 50;
        double c = 0;
        int i = 0;
        int j = 0;
        int k = 0;
        int l = 0;

        double u[][] = Numeric.cloneArray(AtA);
        int m = u.length;

        int n = u[0].length;

        if (m < n) {
            System.err.println("Need more rows than columns");
            return;
        }

        double e[] = new double[n];
        double q[] = new double[n];
        for (i = 0; i < n; i++) {
            e[i] = q[i] = 0.0;
        }
        double v[][] = {{0, 0}, {0, 0}};

        //Householder's reduction to bidiagonal form
        double f = 0.0;
        double g = 0.0;
        double h = 0.0;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        double s = 0.0;

        for (i = 0; i < n; i++) {
            e[i] = g;
            s = 0.0;
            l = i + 1;
            for (j = i; j < m; j++) {
                s += (u[j][i] * u[j][i]);
            }
            if (s <= tolerance) {
                g = 0.0;
            } else {
                f = u[i][i];
                g = Math.sqrt(s);
                if (f >= 0.0) {
                    g = -g;
                }
                h = f * g - s;
                u[i][i] = f - g;
                for (j = l; j < n; j++) {
                    s = 0.0;
                    for (k = i; k < m; k++) {
                        s += u[k][i] * u[k][j];
                    }
                    f = s / h;
                    for (k = i; k < m; k++) {
                        u[k][j] += f * u[k][i];
                    }
                }
            }
            q[i] = g;
            s = 0.0;
            for (j = l; j < n; j++) {
                s = s + u[i][j] * u[i][j];
            }
            if (s <= tolerance) {
                g = 0.0;
            } else {
                f = u[i][i + 1];
                g = Math.sqrt(s);
                if (f >= 0.0) {
                    g = -g;
                }
                h = f * g - s;
                u[i][i + 1] = f - g;
                for (j = l; j < n; j++) {
                    e[j] = u[i][j] / h;
                }
                for (j = l; j < m; j++) {
                    s = 0.0;
                    for (k = l; k < n; k++) {
                        s += (u[j][k] * u[i][k]);
                    }
                    for (k = l; k < n; k++) {
                        u[j][k] += s * e[k];
                    }
                }
            }
            y = Math.abs(q[i]) + Math.abs(e[i]);
            if (y > x) {
                x = y;
            }
        }

        // accumulation of right hand gtransformations
        for (i = n - 1; i != -1; i += -1) {
            if (g != 0.0) {
                h = g * u[i][i + 1];
                for (j = l; j < n; j++) {
                    v[j][i] = u[i][j] / h;
                }
                for (j = l; j < n; j++) {
                    s = 0.0;
                    for (k = l; k < n; k++) {
                        s += u[i][k] * v[k][j];
                    }
                    for (k = l; k < n; k++) {
                        v[k][j] += (s * v[k][i]);
                    }
                }
            }
            for (j = l; j < n; j++) {
                v[i][j] = 0;
                v[j][i] = 0;
            }
            v[i][i] = 1;
            g = e[i];
            l = i;
        }

        // accumulation of left hand transformations
        for (i = n - 1; i != -1; i += -1) {
            l = i + 1;
            g = q[i];
            for (j = l; j < n; j++) {
                u[i][j] = 0;
            }
            if (g != 0.0) {
                h = u[i][i] * g;
                for (j = l; j < n; j++) {
                    s = 0.0;
                    for (k = l; k < m; k++) {
                        s += u[k][i] * u[k][j];
                    }
                    f = s / h;
                    for (k = i; k < m; k++) {
                        u[k][j] += f * u[k][i];
                    }
                }
                for (j = i; j < m; j++) {
                    u[j][i] = u[j][i] / g;
                }
            } else {
                for (j = i; j < m; j++) {
                    u[j][i] = 0;
                }
            }
            u[i][i] += 1;
        }

        // diagonalization of the bidiagonal form
        prec = prec * x;
        for (k = n - 1; k != -1; k += -1) {
            for (int iteration = 0; iteration < itmax; iteration++) {	// test f splitting
                boolean test_convergence = false;
                for (l = k; l != -1; l += -1) {
                    if (Math.abs(e[l]) <= prec) {
                        test_convergence = true;
                        break;
                    }
                    if (Math.abs(q[l - 1]) <= prec) {
                        break;
                    }
                }
                if (!test_convergence) {	// cancellation of e[l] if l>0
                    c = 0.0;
                    s = 1.0;
                    int l1 = l - 1;
                    for (i = l; i < k + 1; i++) {
                        f = s * e[i];
                        e[i] = c * e[i];
                        if (Math.abs(f) <= prec) {
                            break;
                        }
                        g = q[i];
                        h = pythag(f, g);
                        q[i] = h;
                        c = g / h;
                        s = -f / h;
                        for (j = 0; j < m; j++) {
                            y = u[j][l1];
                            z = u[j][i];
                            u[j][l1] = y * c + (z * s);
                            u[j][i] = -y * s + (z * c);
                        }
                    }
                }
                // test f convergence
                z = q[k];
                if (l == k) {	//convergence
                    if (z < 0.0) {	//q[k] is made non-negative
                        q[k] = -z;
                        for (j = 0; j < n; j++) {
                            v[j][k] = -v[j][k];
                        }
                    }
                    break;  //break out of iteration loop and move on to next k value
                }
                if (iteration >= itmax - 1) {
                    System.err.println("Error: no convergence.");
                    return;
                }
                // shift from bottom 2x2 minor
                x = q[l];
                y = q[k - 1];
                g = e[k - 1];
                h = e[k];
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                g = pythag(f, 1.0);
                if (f < 0.0) {
                    f = ((x - z) * (x + z) + h * (y / (f - g) - h)) / x;
                } else {
                    f = ((x - z) * (x + z) + h * (y / (f + g) - h)) / x;
                }
                // next QR transformation
                c = 1.0;
                s = 1.0;
                for (i = l + 1; i < k + 1; i++) {
                    g = e[i];
                    y = q[i];
                    h = s * g;
                    g = c * g;
                    z = pythag(f, h);
                    e[i - 1] = z;
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = -x * s + g * c;
                    h = y * s;
                    y = y * c;
                    for (j = 0; j < n; j++) {
                        x = v[j][i - 1];
                        z = v[j][i];
                        v[j][i - 1] = x * c + z * s;
                        v[j][i] = -x * s + z * c;
                    }
                    z = pythag(f, h);
                    q[i - 1] = z;
                    c = f / z;
                    s = h / z;
                    f = c * g + s * y;
                    x = -s * g + c * y;
                    for (j = 0; j < m; j++) {
                        y = u[j][i - 1];
                        z = u[j][i];
                        u[j][i - 1] = y * c + z * s;
                        u[j][i] = -y * s + z * c;
                    }
                }
                e[l] = 0.0;
                e[k] = f;
                q[k] = x;
            }
        }

        //vt= transpose(v)
        //return (u,q,vt)
        for (i = 0; i < q.length; i++) {
            if (q[i] < prec) {
                q[i] = 0;
            }
        }

        //sort eigenvalues	
        for (i = 0; i < n; i++) {
            //writeln(q)
            for (j = i - 1; j >= 0; j--) {
                if (q[j] < q[i]) {
                    //  writeln(i,'-',j)
                    c = q[j];
                    q[j] = q[i];
                    q[i] = c;
                    for (k = 0; k < u.length; k++) {
                        temp = u[k][i];
                        u[k][i] = u[k][j];
                        u[k][j] = temp;
                    }
                    for (k = 0; k < v.length; k++) {
                        temp = v[k][i];
                        v[k][i] = v[k][j];
                        v[k][j] = temp;
                    }
//	   u.swapCols(i,j)
//	   v.swapCols(i,j)
                    i = j;
                }
            }
        }

        //return {U:u,S:q,V:v}
        // SVD vége
        double Ut[][] = {{u[0][0], u[1][0]}, {u[0][1], u[1][1]}};
        double D[][] = {{1 / Math.sqrt(q[0]), 0}, {0, 1 / Math.sqrt(q[1])}};
        double Sinv[][] = {
            {
                (D[0][0] * Ut[0][0] + D[0][1] * Ut[1][0]) * v[0][0] + (D[1][0] * Ut[0][0] + D[1][1] * Ut[1][0]) * v[0][1],
                (D[0][0] * Ut[0][1] + D[0][1] * Ut[1][1]) * v[0][0] + (D[1][0] * Ut[0][1] + D[1][1] * Ut[1][1]) * v[0][1]
            },
            {
                (D[0][0] * Ut[0][0] + D[0][1] * Ut[1][0]) * v[1][0] + (D[1][0] * Ut[0][0] + D[1][1] * Ut[1][0]) * v[1][1],
                (D[0][0] * Ut[0][1] + D[0][1] * Ut[1][1]) * v[1][0] + (D[1][0] * Ut[0][1] + D[1][1] * Ut[1][1]) * v[1][1]
            }
        };
        double R[][] = {
            {A[0][0] * Sinv[0][0] + A[0][1] * Sinv[1][0], A[0][0] * Sinv[0][1] + A[0][1] * Sinv[1][1]},
            {A[1][0] * Sinv[0][0] + A[1][1] * Sinv[1][0], A[1][0] * Sinv[0][1] + A[1][1] * Sinv[1][1]}
        };

        double stiffness_adapted = 1 - Math.pow(1 - this.stiffness, 1 / Physics.solver_iterations);

        for (ParticleSMC p : particles) {
            double goal_position[] = {
                com[0] + R[0][0] * p.rxRest + R[0][1] * p.ryRest,
                com[1] + R[1][0] * p.rxRest + R[1][1] * p.ryRest
            };
            double delta_p[] = {
                goal_position[0] - p.tempX,
                goal_position[1] - p.tempY
            };
            p.tempX += stiffness_adapted * delta_p[0];
            p.tempY += stiffness_adapted * delta_p[1];
        }

        int a = 0;
    }

    private double pythag(double pa, double pb) {

        double a = Math.abs(pa);
        double b = Math.abs(pb);

        if (a > b) {
            return a * Math.sqrt(1.0 + (b * b / a / a));
        } else if (b == 0.0) {
            return a;
        }

        return b * Math.sqrt(1.0 + (a * a / b / b));
    }

}
