package workshop;

import Jama.Matrix;
import jv.geom.PgElementSet;
import jv.object.PsDebug;
import jv.vecmath.PdMatrix;
import jv.vecmath.PdVector;
import util.Util;

import java.util.Arrays;
import java.util.Collections;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class RigidRegistration {
    public PgElementSet p;
    public PgElementSet q;
    private int k;
    private boolean pointToPlane;

    PdMatrix optRotation;
    PdMatrix optTranslation;

    public RigidRegistration(PgElementSet p, PgElementSet q, int k, boolean pointToPlane) {
        this.p = p;
        this.q = q;
        this.k = k;
        this.pointToPlane = pointToPlane;
    }

    // get average point from a list of vectors
    private PdVector average(PdVector[] v) {
        PdVector sum = Arrays.stream(v).reduce(PdVector::addNew).get();
        sum.multScalar(1.0 / v.length);
        return sum;

    }

    public void runAlgorithm() {
        if (pointToPlane) {
            pointToPlane();
        } else {
            pointToPoint();
        }
    }

    private void pointToPlane() {
        // step 0: closest points
        Set<VertexPair> closestPairs = trimPairs(computeClosestPairs(this::distanceToPlane), k);

        // compute a & b
        PdMatrix A = new PdMatrix(6, 6);
        PdMatrix b = new PdMatrix(6, 1);

        PsDebug.message("Running loop: " + closestPairs.size());
        for (VertexPair pair : closestPairs) {
            PdVector normal = q.getVertexNormal(pair.v1);
            PdVector pointP = p.getVertex(pair.v0);

            // step 1: compute A
            PdMatrix topRow = toMatrix(PdVector.crossNew(pointP, normal));
            PdMatrix bottomRow = toMatrix(normal);

            PdMatrix ai = new PdMatrix(6, 1);
            ai.set(new double[]{
                    topRow.getRow(0).getEntry(0),
                    topRow.getRow(1).getEntry(0),
                    topRow.getRow(2).getEntry(0),
                    bottomRow.getRow(0).getEntry(0),
                    bottomRow.getRow(1).getEntry(0),
                    bottomRow.getRow(2).getEntry(0)
            });

            PdMatrix ai_t = new PdMatrix(1, 6);
            ai_t.transpose(ai);

            PdMatrix res = new PdMatrix(6, 6);
            res.mult(ai, ai_t);

            A.add(res);

            // step 2: compute B
            double distance = signedPlaneDistance(pair.v0, pair.v1);

            PdMatrix resB = new PdMatrix(6, 1);
            resB.multScalar(ai, distance);

            b.add(resB);
        }

        // step 3: solve ...
        PdMatrix minusB = new PdMatrix(6, 1);
        minusB.set(b.m_data);
        minusB.multScalar(-1);

        Matrix jama = new Matrix(A.m_data);
        Matrix rt = jama.solve(new Matrix(minusB.m_data));
        PdMatrix r = new PdMatrix(3, 1);
        r.set(new double[]{
                rt.get(0, 0),
                rt.get(1, 0),
                rt.get(2, 0)
        });
        PdMatrix t = new PdMatrix(3, 1);
        t.set(new double[]{
                rt.get(3, 0),
                rt.get(4, 0),
                rt.get(5, 0)
        });

        double r1 = r.getRow(0).getEntry(0);
        double r2 = r.getRow(1).getEntry(0);
        double r3 = r.getRow(2).getEntry(0);
        PdMatrix R_apos = new PdMatrix(new double[][]{
                {1, -r3, r2},
                {r3, 1, -r1},
                {-r2, r1, 1}
        });

        PdMatrix u = new PdMatrix(3, 3);
        PdMatrix d = new PdMatrix(3, 3);
        PdMatrix v = new PdMatrix(3, 3);
        Util.computeSVD(R_apos, u, d, v);

        PdMatrix vt = transpose(v);

        // step 4: compute optional rotation
        PsDebug.message("Running step 4");
        PdMatrix rOpt = new PdMatrix(3, 3);

        PdMatrix uvt = new PdMatrix(3, 3);
        uvt.mult(u, vt);

        PdMatrix middle = new PdMatrix(new double[][]{
                {1.0, 0, 0},
                {0, 1, 0},
                {0, 0, uvt.det()}
        });

        PdMatrix rOptIntermediateResult = new PdMatrix(3, 3);
        rOptIntermediateResult.mult(middle, vt);
        rOpt.mult(u, rOptIntermediateResult);
        PsDebug.message("ropt:" + rOpt.toString());
        PsDebug.message("topt:" + t.toString());

        transform(toMatrix(t.getColumn(0)), rOpt);
    }

    private void pointToPoint() {
        // step 0: closest points
        Set<VertexPair> closestPairs = trimPairs(computeClosestPairs(this::distance), k);

        // step 1: compute centroids
        PsDebug.message("Running step 1");
        PdVector centroidP = average(p.getVertices());
        PdVector centroidQ = average(q.getVertices());

        // step 2: compute M = ...
        PsDebug.message("Running step 2");
        PdMatrix m = new PdMatrix(3, 3);
        int n = closestPairs.size();

        for (VertexPair pair : closestPairs) {
            PdMatrix dp = new PdMatrix(3, 1);
            PdMatrix dq = new PdMatrix(3, 1);
            PdMatrix dqt = new PdMatrix(1, 3);

            dp.set(PdVector.subNew(p.getVertex(pair.v0), centroidP).m_data);
            dq.set(PdVector.subNew(q.getVertex(pair.v1), centroidQ).m_data);

            dqt.transpose(dq);

            PdMatrix res = new PdMatrix(3, 3);
            res.mult(dp, dqt);
            res.multScalar(1.0 / n);

            m.add(res);
        }

        // step 3: compute SVD
        PsDebug.message("Running step 3");
        PdMatrix u = new PdMatrix(3, 3);
        PdMatrix d = new PdMatrix(3, 3);
        PdMatrix v = new PdMatrix(3, 3);

        Util.computeSVD(m, u, d, v);

        PdMatrix ut = new PdMatrix(3, 3);
        ut.transpose(u);

        // step 4: compute optional rotation
        PsDebug.message("Running step 4");
        PdMatrix rOpt = new PdMatrix(3, 3);

        PdMatrix vut = new PdMatrix(3, 3);
        vut.mult(v, ut);

        PdMatrix middle = new PdMatrix(new double[][]{
                {1.0, 0, 0},
                {0, 1, 0},
                {0, 0, vut.det()}
        });

        PdMatrix rOptIntermediateResult = new PdMatrix(3, 3);
        rOptIntermediateResult.mult(middle, ut);
        rOpt.mult(v, rOptIntermediateResult);

        optRotation = rOpt;

        // step 5: compute optimal translation
        PsDebug.message("Running step 5");
        PdMatrix centroidPMatrix = new PdMatrix(3, 1);
        PdMatrix centroidQMatrix = new PdMatrix(3, 1);
        centroidPMatrix.set(centroidP.m_data);
        centroidQMatrix.set(centroidQ.m_data);

        PsDebug.message(centroidPMatrix.toString());
        PsDebug.message(rOpt.toString());

        PdMatrix rOptIntermediate = new PdMatrix(3, 1);
        rOptIntermediate.mult(rOpt, centroidPMatrix);
        centroidQMatrix.sub(rOptIntermediate);

        optTranslation = centroidQMatrix;

        PsDebug.message(this.toString());

        transform(optTranslation, optRotation);
    }


    private void transform(PdMatrix translation, PdMatrix rotation) {
        PdMatrix res = new PdMatrix(3, 1);
        for (PdVector v : p.getVertices()) {
            res.mult(rotation, toMatrix(v));
            res.add(translation);

            v.set(res.getColumn(0).m_data);
        }

        p.update(p);
    }

    @Override
    public String toString() {
        return "RigidRegistration{" +
                "optRotation=" + optRotation +
                ", optTranslation=" + optTranslation +
                '}';
    }

    private Set<VertexPair> trimPairs(Set<VertexPair> pairs, int k) {
        // https://stackoverflow.com/a/49215170
        double median = pairs.stream().mapToDouble(VertexPair::getDistance).sorted()
                .skip((pairs.size() - 1) / 2).limit(2 - pairs.size() % 2).average().getAsDouble();

        PsDebug.message("Median distance: " + median);

        return pairs.stream().filter(p -> p.getDistance() <= k * median).collect(Collectors.toSet());
    }

    // since finding closest pairs is slow without a fancy data structure, use a ConcurrentSet and a parallel stream to
    // do it concurrently
    public Set<VertexPair> computeClosestPairs(DistanceFunction fnDistance) {
        Set<VertexPair> res = Collections.newSetFromMap(new ConcurrentHashMap<>());

        IntStream.range(0, p.getNumVertices()).parallel().forEach(i -> {
            res.add(new VertexPair(i, findClosestInQ(i, fnDistance), fnDistance));
        });

        return res;
    }

    private int findClosestInQ(int v, DistanceFunction fnDistance) {
        int pMin = -1;
        double distMin = Double.MAX_VALUE;

        PdVector[] vertices = q.getVertices();
        for (int i = 0, verticesLength = vertices.length; i < verticesLength; i++) {
            // since we're just comparing, square distance is enough
            double dist = fnDistance.apply(v, i);
            if (dist < distMin) {
                distMin = dist;
                pMin = i;
            }
        }
        return pMin;
    }

    private double squareDistance(int v0, int v1) {
        PdVector diff = PdVector.subNew(q.getVertex(v0), p.getVertex(v1));
        return PdVector.dot(diff, diff);
    }

    private double signedPlaneDistance(int v0, int v1) {
        PdVector pointP = p.getVertex(v0);
        PdVector pointQ = q.getVertex(v1);
        PdVector normalQ = q.getVertexNormal(v1);

        PdMatrix difference = toMatrix(PdVector.subNew(pointP, pointQ));
        PdMatrix differenceTranspose = transpose(difference);

        PdMatrix res = new PdMatrix(1, 1);

        res.mult(differenceTranspose, toMatrix(normalQ));

        return res.getEntry(0, 0);
    }

    private double distanceToPlane(int v0, int v1) {
        return Math.abs(signedPlaneDistance(v0, v1));
    }

    private double distance(int v0, int v1) {
        return Math.sqrt(squareDistance(v0, v1));
    }


    class VertexPair {
        int v0;
        int v1;

        double distance;

        public VertexPair(int v0, int v1, DistanceFunction distFn) {
            this.v0 = v0;
            this.v1 = v1;

            distance = distFn.apply(v0, v1);
        }


        public double getDistance() {
            return distance;
        }

        @Override
        public String toString() {
            return "VertexPair{" +
                    "v0=" + Arrays.toString(p.getVertex(v0).m_data) +
                    ", v1=" + Arrays.toString(q.getVertex(v1).m_data) +
                    ", distance=" + distance +
                    '}';
        }
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
}

interface DistanceFunction {
    double apply(int a, int b);
}

