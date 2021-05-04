package workshop;

import Jama.Matrix;
import jv.geom.PgElementSet;
import jv.object.PsDebug;
import jv.vecmath.PdMatrix;
import jv.vecmath.PdMatrixIf;
import jv.vecmath.PdVector;
import jv.vecmath.PuMath;
import util.Util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

public class RigidRegistration {
    public PgElementSet p;
    public PgElementSet q;

    Set<VertexPair> closestPairs;

    PdMatrix optRotation;
    PdMatrix optTranslation;

    public RigidRegistration(PgElementSet p, PgElementSet q, int k) {
        this.p = p;
        this.q = q;

        // part 1
        Set<VertexPair> pairs = computeClosestPairs();

        // part 2
        closestPairs = trimPairs(pairs, k);
    }

    // get average point from a list of vectors
    private PdVector average(PdVector[] v) {
        PdVector sum = Arrays.stream(v).reduce(PdVector::addNew).get();
        sum.multScalar(1.0 / v.length);
        return sum;
    }

    public void runAlgorithm() {
        // step 1: compute centroids
        PsDebug.message("Running step 1");
        PdVector centroidP = average(p.getVertices());
        PdVector centroidQ = average(q.getVertices());

        // step 2: compute M = ...
        PdMatrix m = new PdMatrix(3, 3);
        int n = closestPairs.size();

        PsDebug.message("Running step 2");
        for (VertexPair pair : closestPairs) {
            PdMatrix dp = new PdMatrix(3, 1);
            PdMatrix dq = new PdMatrix(3, 1);

            dp.set(PdVector.subNew(pair.v0, centroidP).m_data);
            dq.set(PdVector.subNew(pair.v1, centroidQ).m_data);
            dq.transpose();

            PdMatrix res = new PdMatrix(3, 3);
            res.mult(dp, dq);
            res.multScalar(1.0 / n);
            m.add(res);
        }

        // step 3: compute SVD
        PsDebug.message("Running step 3");
        PdMatrix u = new PdMatrix(3, 3);
        PdMatrix d = new PdMatrix(3, 3);
        PdMatrix v = new PdMatrix(3, 3);

        PdMatrix ut = new PdMatrix(3, 3);
        ut.transpose(u);

        Util.computeSVD(m, u, d, v);

        // step 4: compute optional rotation
        PsDebug.message("Running step 4");
        PdMatrix rOpt = new PdMatrix(3, 3);
        PdMatrix rOptIntermediateResult = new PdMatrix(3, 3);

        PdMatrix vut = new PdMatrix(3, 3);
        vut.mult(v, ut);

        PdMatrix middle = new PdMatrix(new double[][]{
                {1.0, 0, 0},
                {0, 1, 0},
                {0, 0, vut.det()}
        });

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
        q.translate(translation.getColumn(0));

        for (PdVector v : q.getVertices()) {
            v.leftMultMatrix(rotation);
        }
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

        return pairs.stream().filter(p -> p.getDistance() > k * median).collect(Collectors.toSet());
    }

    // since finding closest pairs is slow without a fancy data structure, use a ConcurrentSet and a parallel stream to
    // do it concurrently
    public Set<VertexPair> computeClosestPairs() {
        Set<VertexPair> res = Collections.newSetFromMap(new ConcurrentHashMap<>());

        Arrays.stream(q.getVertices()).parallel().forEach(v -> {
            res.add(new VertexPair(v, findClosestInQ(v)));
        });

        return res;
    }

    private PdVector findClosestInQ(PdVector p0) {
        PdVector pMin = null;
        double distMin = Double.MAX_VALUE;

        for (PdVector p1 : q.getVertices()) {
            // since we're just comparing, square distance is enough
            double dist = squareDistance(p0, p1);
            if (dist < distMin) {
                distMin = dist;
                pMin = p1;
            }
        }
        return pMin;
    }

    private double squareDistance(PdVector v0, PdVector v1) {
        PdVector diff = PdVector.subNew(v0, v1);
        return PdVector.dot(diff, diff);
    }

    private double distance(PdVector v0, PdVector v1) {
        return Math.sqrt(squareDistance(v0, v1));
    }


    class VertexPair {
        PdVector v0;
        PdVector v1;

        double distance;

        public VertexPair(PdVector v0, PdVector v1) {
            this.v0 = v0;
            this.v1 = v1;

            distance = distance(v0, v1);
        }

        public double getDistance() { return distance; }
    }
}

