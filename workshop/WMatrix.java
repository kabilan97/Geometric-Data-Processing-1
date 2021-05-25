package workshop;

import Jama.Matrix;
import jv.vecmath.PdMatrix;

import java.util.Arrays;

public class WMatrix {
    PdMatrix m;

    public WMatrix(double[][] v) {
        this.m = new PdMatrix(v);
    }

    public WMatrix(PdMatrix m) {
        this.m = m;
    }

    public WMatrix(int rows, int cols) {
        this.m = new PdMatrix(rows, cols);
    }

    public static WMatrix columns(WVector... col) {
        PdMatrix m = new PdMatrix(col[0].v.getSize(), col.length);
        for (int i = 0; i < col.length; i++) {
            m.setColumn(i, col[i].v);
        }

        return new WMatrix(m);
    }

    public WMatrix mult(WMatrix o) {
        assert m.getNumCols() == o.m.getNumRows() : "Matrix dimensions don't allow for multiplication";

        PdMatrix res = new PdMatrix(m.getNumRows(), o.m.getNumCols());
        res.mult(m, o.m);

        return new WMatrix(res);
    }

    @Override
    public String toString() {
        return Arrays.deepToString(m.m_data).replace("], ", "]\n");
    }

    public WMatrix mult(WVector n) {
        return mult(n.toMatrix());
    }

    public double asNumber() {
        assert m.getNumRows() == 1 && m.getNumCols() == 1 : "Matrix is not 1x1: " + m.getNumRows() + "x" + m.getNumCols();

        return m.getEntry(0, 0);
    }

    public WMatrix add(WMatrix o) {
        PdMatrix res = new PdMatrix(m.getNumRows(), m.getNumCols());
        res.add(this.m, o.m);
        return new WMatrix(res);
    }


    public WMatrix transpose() {
        PdMatrix res = new PdMatrix(m.getNumCols(), m.getNumRows());
        res.transpose(m);
        return new WMatrix(res);
    }

    public WMatrix solve(WMatrix other) {
        Matrix jama = new Matrix(m.m_data);
        Matrix solved = jama.solve(new Matrix(other.m.m_data));

        return new WMatrix(solved.getArray());
    }

    public WMatrix solve(WVector vector) {
        return solve(vector.toMatrix());
    }

    public WMatrix mult(double scalar) {
        PdMatrix res = new PdMatrix(m.getNumRows(), m.getNumCols());
        res.multScalar(m, scalar);
        return new WMatrix(res);
    }

    public WMatrix negative() {
        return mult(-1);
    }

    public WVector toVector() {
        return new WVector(m.getColumn(0));
    }
}
