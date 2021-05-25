package workshop;

import jv.object.PsDebug;
import jv.vecmath.PdMatrix;
import jv.vecmath.PdVector;

import java.util.Arrays;

public class WVector {
    PdVector v;

    public WVector(double... a) {
        this.v = new PdVector(a);
    }

    public WVector(PdVector v) {
        this.v = v;
    }

    public WMatrix transpose() {
        return new WMatrix(transpose(toMatrix().m));
    }

    public double x() {
        return v.getEntry(0);
    }

    public double y() {
        return v.getEntry(1);
    }

    public double z() {
        return v.getEntry(2);
    }

    public WMatrix toMatrix() {
        return new WMatrix(toMatrix(v));
    }

    public WMatrix mult(WMatrix other) {
        return this.toMatrix().mult(other);
    }

    public double dot(WVector o) {
        return v.dot(o.v);
    }

    public WVector minus(WVector o) {
        PdVector res = new PdVector();
        res.sub(this.v, o.v);

        return new WVector(res);
    }

    public WVector plus(WVector o) {
        PdVector res = new PdVector();
        res.add(this.v, o.v);

        return new WVector(res);
    }


    private static PdMatrix toMatrix(PdVector v) {
        PdMatrix m = new PdMatrix(v.getSize(), 1);
        m.set(v.m_data);

        return m;
    }

    private static PdMatrix transpose(PdMatrix m) {
        PdMatrix r = new PdMatrix(m.getNumCols(), m.getNumRows());
        r.transpose(m);
        return r;
    }

    @Override
    public String toString() {
        return Arrays.toString(v.m_data);
    }

    public WVector cross(WVector other) {
        return new WVector(PdVector.crossNew(this.v, other.v));
    }

    public double length() {
        return v.length();
    }
}
