package workshop;

import Jama.Matrix;
import jv.geom.PgElementSet;
import jv.object.PsDebug;
import jv.vecmath.PdMatrix;
import jv.vecmath.PdVector;
import util.Util;

import java.util.Arrays;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class RigidRegistration {
    public PgElementSet p;
    public PgElementSet q;
    private int k;
    private int n;
    private boolean pointToPlane;

    double totalTranslation = 0.0;

    public double getTotalTranslation() {
        return totalTranslation;
    }

    public RigidRegistration(PgElementSet p, PgElementSet q, int n, int k, boolean pointToPlane) {
        this.p = p;
        this.q = q;
        this.n = n;
        this.k = k;
        this.pointToPlane = pointToPlane;
    }

    // get average point from a list of vectors
    private PdVector average(PdVector[] v) {
        PdVector sum = Arrays.stream(v).reduce(PdVector::addNew).get();
        sum.multScalar(1.0 / v.length);
        return sum;
    }

    private PdVector average(Set<PdVector> v) {
        PdVector sum = v.stream().reduce(PdVector::addNew).get();
        sum.multScalar(1.0 / v.size());
        return sum;
    }


    public void runAlgorithm() {
        try {
            if (pointToPlane) {
                PsDebug.message("Point to Plane");
                pointToPlane();
            } else {
                PsDebug.message("Point to Point");
                pointToPoint();
            }
        } catch (Exception ex) {
            PsDebug.message(Arrays.toString(ex.getStackTrace()));
        }
    }

    private void pointToPlane() {
        // step 0: closest points
        Set<VertexPair> closestPairs = computeClosestPairs(this::squareDistance);

        // compute a & b
        PdMatrix A = new PdMatrix(6, 6);
        PdMatrix b = new PdMatrix(6, 1);

        PsDebug.message("Running loop: " + closestPairs.size());
        for (VertexPair pair : closestPairs) {
            PdVector normal = q.getVertexNormal(pair.indexQ);
            PdVector pointP = p.getVertex(pair.indexP);

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
            double distance = signedPlaneDistance(pair.indexP, pair.indexQ);

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
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, uvt.det()}
        });

        PdMatrix rOptIntermediateResult = new PdMatrix(3, 3);
        rOptIntermediateResult.mult(u, middle);
        rOpt.mult(rOptIntermediateResult, vt);

        transform(toMatrix(t.getColumn(0)), rOpt);
    }

    private void pointToPoint() {
        // step 0: closest points
        Set<VertexPair> closestPairs = computeClosestPairs(this::squareDistance);

        // step 1: compute centroids
        PsDebug.message("Running step 1");

//        PdVector centroidP = average(p.getVertices());
//        PdVector centroidQ = average(q.getVertices());
        PdVector centroidP = average(closestPairs.stream().map(pair -> p.getVertex(pair.indexP)).collect(Collectors.toSet()));
        PdVector centroidQ = average(closestPairs.stream().map(pair -> q.getVertex(pair.indexQ)).collect(Collectors.toSet()));

//        PsDebug.message("Old: " + Arrays.toString(centroidPOld.m_data));
//        PsDebug.message("New: " + Arrays.toString(centroidP.m_data));

        // step 2: compute M = ...
        PsDebug.message("Running step 2");
        PdMatrix m = new PdMatrix(3, 3);
        int n = closestPairs.size();

        for (VertexPair pair : closestPairs) {
            PdMatrix dp = new PdMatrix(3, 1);
            PdMatrix dq = new PdMatrix(3, 1);
            PdMatrix dqt = new PdMatrix(1, 3);

            dp.set(PdVector.subNew(p.getVertex(pair.indexP), centroidP).m_data);
            dq.set(PdVector.subNew(q.getVertex(pair.indexQ), centroidQ).m_data);

            dqt.transpose(dq);

            PdMatrix res = new PdMatrix(3, 3);
            res.mult(dp, dqt);

            m.add(res);
        }
        m.multScalar(1.0 / n);

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
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, vut.det()}
        });

        PdMatrix rOptIntermediateResult = new PdMatrix(3, 3);
        rOptIntermediateResult.mult(v, middle);
        rOpt.mult(rOptIntermediateResult, ut);

        PdMatrix optRotation = rOpt;

        // step 5: compute optimal translation
        PsDebug.message("Running step 5");
        PdMatrix centroidPMatrix = new PdMatrix(3, 1);
        PdMatrix centroidQMatrix = new PdMatrix(3, 1);
        centroidPMatrix.set(centroidP.m_data);
        centroidQMatrix.set(centroidQ.m_data);

        PdMatrix rOptIntermediate = new PdMatrix(3, 1);
        rOptIntermediate.mult(rOpt, centroidPMatrix);

        PdMatrix optTranslation = new PdMatrix(3, 1);
        optTranslation.sub(centroidQMatrix, rOptIntermediate);

        transform(optTranslation, optRotation);
    }


    private void transform(PdMatrix translation, PdMatrix rotation) {
        PsDebug.message("ropt:" + rotation.toString());
        PsDebug.message("topt:" + translation.toString());

        PdMatrix res = new PdMatrix(3, 1);
        for (PdVector v : p.getVertices()) {
            res.mult(rotation, toMatrix(v));
            res.add(translation);

            v.set(res.getColumn(0).m_data);
        }

        PdVector translationVector = new PdVector(translation.getColumn(0).m_data);
        this.totalTranslation = translationVector.length();

        p.update(p);
    }

    // since finding closest pairs is slow without a fancy data structure, parallel stream to
    // do it concurrently
    public Set<VertexPair> computeClosestPairs(DistanceFunction fnDistance) {
        IntStream stream;
        if (n < 1 || n > p.getNumVertices()) {
            stream = IntStream.range(0, p.getNumVertices());
        } else {
            stream = new Random().ints(n,0, p.getNumVertices()).distinct();
        }

        Set<VertexPair> pairs = stream.parallel()
                .mapToObj(indexP -> new VertexPair(indexP, findClosestInQ(indexP, fnDistance), fnDistance))
                .collect(Collectors.toSet());

        // https://stackoverflow.com/a/49215170
        double median = pairs.stream().parallel().mapToDouble(VertexPair::getDistance).sorted()
                .skip((pairs.size() - 1) / 2).limit(2 - pairs.size() % 2).average().getAsDouble();

        PsDebug.message("Median distance: " + median);

        Set<VertexPair> trimmed = pairs.stream().filter(p -> p.getDistance() <= k * median)
                .collect(Collectors.toSet());

        return trimmed;
    }

    private int findClosestInQ(int indexP, DistanceFunction fnDistance) {
        int minPoint = -1;
        double minDist = Double.MAX_VALUE;

        PdVector[] vertices = q.getVertices();
        for (int indexQ = 0, verticesLength = vertices.length; indexQ < verticesLength; indexQ++) {
            // since we're just comparing, square distance is enough
            double dist = fnDistance.apply(indexP, indexQ);
            if (dist < minDist) {
                minDist = dist;
                minPoint = indexQ;
            }
        }
        return minPoint;
    }

    private double squareDistance(int v0, int v1) {
        PdVector diff = PdVector.subNew(p.getVertex(v0), q.getVertex(v1));
        return PdVector.dot(diff, diff);
    }

    private double signedPlaneDistance(int indexP, int indexQ) {
        PdVector pointP = p.getVertex(indexP);
        PdVector pointQ = q.getVertex(indexQ);
        PdVector normalQ = q.getVertexNormal(indexQ);

        PdMatrix difference = toMatrix(PdVector.subNew(pointP, pointQ));
        PdMatrix differenceTranspose = transpose(difference);

        PdMatrix res = new PdMatrix(1, 1);
        res.mult(differenceTranspose, toMatrix(normalQ));

        return res.getEntry(0, 0);
    }

    private double distanceToPlane(int indexP, int indexQ) {
        return Math.sqrt(signedPlaneDistance(indexP, indexQ));
    }

    private double distance(int v0, int v1) {
        return Math.sqrt(squareDistance(v0, v1));
    }


    class VertexPair {
        int indexP;
        int indexQ;

        double distance;

        public VertexPair(int indexP, int indexQ, DistanceFunction distFn) {
            this.indexP = indexP;
            this.indexQ = indexQ;

            distance = distFn.apply(indexP, indexQ);
        }


        public double getDistance() {
            return distance;
        }

        @Override
        public String toString() {
            return "VertexPair{" +
                    "v0=" + Arrays.toString(p.getVertex(indexP).m_data) +
                    ", v1=" + Arrays.toString(q.getVertex(indexQ).m_data) +
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

