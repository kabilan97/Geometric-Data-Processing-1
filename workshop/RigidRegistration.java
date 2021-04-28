package workshop;

import jv.geom.PgElementSet;
import jv.vecmath.PdVector;

import java.util.Arrays;
import java.util.Collections;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

public class RigidRegistration {
    public PgElementSet p;
    public PgElementSet q;

    public RigidRegistration(PgElementSet p, PgElementSet q, int k) {
        this.p = p;
        this.q = q;

        // part 1
        Set<VertexPair> pairs = computeClosestPairs();

        // part 2
        Set<VertexPair> trimmedPairs = trimPairs(pairs, k);

        // part 3
        // TODO
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

