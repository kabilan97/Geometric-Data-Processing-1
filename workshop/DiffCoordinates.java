package workshop;

import jv.geom.PgElementSet;
import jv.object.PsDebug;
import jv.object.PsObject;
import jv.project.PgGeometry;
import jv.vecmath.PdMatrix;
import jv.vecmath.PdVector;
import jv.vecmath.PiVector;
import jvx.numeric.PnConjugateGradient;
import jvx.numeric.PnConjugateGradientMatrix;
import jvx.numeric.PnSparseMatrix;
import jvx.project.PjWorkshop;
import dev6.numeric.PnMumpsSolver;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import static jvx.numeric.PnSparseMatrix.*;

public class DiffCoordinates extends PjWorkshop {

	PgElementSet mesh;
	PgElementSet m_geomSave;

	public DiffCoordinates() {
		super("Differential Coordinates");
		init();
	}
	
	@Override
	public void setGeometry(PgGeometry geom) {
		super.setGeometry(geom);
		mesh = (PgElementSet)super.m_geom;
		m_geomSave 	= (PgElementSet)super.m_geomSave;
	}
	
	public void init() {		
		super.init();
	}


	private Matrices computeMatrices() {
		PnSparseMatrix G = new PnSparseMatrix(3 * mesh.getNumElements(), mesh.getNumVertices());
		PnSparseMatrix Mv = new PnSparseMatrix(3 * mesh.getNumElements(), 3 * mesh.getNumElements());
		PnSparseMatrix M = new PnSparseMatrix(mesh.getNumVertices(), mesh.getNumVertices());

		PiVector[] elements = mesh.getElements();
		for (int triangleIndex = 0, elementsLength = elements.length; triangleIndex < elementsLength; triangleIndex++) {
			// Get the triangle vertices
			PiVector vertexIndices = elements[triangleIndex];
			WVector p1 = new WVector(mesh.getVertex(vertexIndices.getEntry(0)));
			WVector p2 = new WVector(mesh.getVertex(vertexIndices.getEntry(1)));
			WVector p3 = new WVector(mesh.getVertex(vertexIndices.getEntry(2)));

			// Compute matrix G
			WVector N = new WVector(mesh.getElementNormal(triangleIndex));
			WMatrix A = computeA(N, p1, p2, p3);
			for (int col = 0; col < 3; col++) {
				for (int row = 0; row < 3; row++) {
					G.setEntry(triangleIndex * 3 + row, vertexIndices.getEntry(col), A.m.getEntry(row, col));
				}
			}

			// Compute matrix Mv
			double area = area(p1, p2, p3);
			for (int i = 0; i < 3; i++) {
				Mv.setEntry(triangleIndex * 3 + i,triangleIndex * 3 + i, area);
			}

			// add area to M
			for (int i = 0; i < 3; i++) {
				int n = vertexIndices.getEntry(i);

				double current = M.getEntry(n, n);
				M.setEntry(n, n, current + area / 3.0);
			}
		}

		PnSparseMatrix Gt = G.transposeNew();
		PnSparseMatrix S = multMatrices(Gt, multMatrices(Mv, G, null), null);
		PnSparseMatrix Mi = new PnSparseMatrix(mesh.getNumVertices(), mesh.getNumVertices());
		for (int i = 0; i < mesh.getNumVertices(); i++) {
			Mi.setEntry(i, i, 1.0 / M.getEntry(i, i));
		}

		PnSparseMatrix L = multMatrices(Mi, S, null);

		return new Matrices(G, Gt, Mv, S, L, L.transposeNew(), M, Mi);
	}

	public void deformMesh(WMatrix userSuppliedMatrix, Set<Integer> selected) throws Exception {
		mesh.makeElementNormals();

		Matrices m = computeMatrices();

		// Step 1: build vector g
		WVector[] v = computeV();
		WVector vx = v[0];
		WVector vy = v[1];
		WVector vz = v[2];

		WVector gx = new WVector(rightMultVector(m.G, vx.v, null));
		WVector gy = new WVector(rightMultVector(m.G, vy.v, null));
		WVector gz = new WVector(rightMultVector(m.G, vz.v, null));

		// Step 2: build gTilde vectors
		WVector gxTilde = new WVector(new PdVector(gx.v.getSize()));
		WVector gyTilde = new WVector(new PdVector(gy.v.getSize()));
		WVector gzTilde = new WVector(new PdVector(gz.v.getSize()));

		for (int i = 0; i < gx.v.getSize(); i++) {
			WVector a = new WVector(new PdVector(gx.v.getEntry(i), gy.v.getEntry(i), gz.v.getEntry(i)));
			WVector res;

			int triangleIndex = i / 3;
			if (selected.contains(mesh.getElement(triangleIndex).getEntry(i % 3))) {
				res = userSuppliedMatrix.mult(a).toVector();
			} else {
				res = a;
			}

			gxTilde.v.setEntry(i, res.x());
			gyTilde.v.setEntry(i, res.y());
			gzTilde.v.setEntry(i, res.z());
		}

		// Step 3: make LHS matrix
		PnSparseMatrix lhs = multMatrices(m.Gt, multMatrices(m.Mv, m.G, null), null);
		lhs.validate();

		// Step 4: make RHS matrix
		PdVector rhsx = rightMultVector(m.Gt, rightMultVector(m.Mv, gxTilde.v, null), null);
		PdVector rhsy = rightMultVector(m.Gt, rightMultVector(m.Mv, gyTilde.v, null), null);
		PdVector rhsz = rightMultVector(m.Gt, rightMultVector(m.Mv, gzTilde.v, null), null);

		// Step 5: Solve to get vector v
		applyAfterSolving(lhs, rhsx, rhsy, rhsz);
	}

	public void deformMeshLaplacian(WMatrix userSuppliedMatrix, Set<Integer> selected, Set<Integer> constrained) throws Exception {
		mesh.makeElementNormals();
		double lambda = 1.0;

		Matrices m = computeMatrices();

		// Step 1: build vectors delta
		WVector[] v = computeV();
		WVector vx = v[0];
		WVector vy = v[1];
		WVector vz = v[2];

		// 1 0 0
		// 0 -0.5 -0.86
		// 0 0.86 -0.5


		WVector deltagx = new WVector(rightMultVector(m.L, vx.v, null));
		WVector deltagy = new WVector(rightMultVector(m.L, vy.v, null));
		WVector deltagz = new WVector(rightMultVector(m.L, vz.v, null));

		PdVector deltax = new PdVector(mesh.getNumVertices());
		PdVector deltay = new PdVector(mesh.getNumVertices());
		PdVector deltaz = new PdVector(mesh.getNumVertices());

		for (int i = 0; i < deltagx.v.getSize(); i++) {
			WVector a = new WVector(new PdVector(deltagx.v.getEntry(i), deltagy.v.getEntry(i), deltagz.v.getEntry(i)));
			WVector res;

			if (selected.contains(i)) {
				res = a;
			} else {
				res = userSuppliedMatrix.mult(a).toVector();
			}

			deltax.setEntry(i, res.x());
			deltay.setEntry(i, res.y());
			deltaz.setEntry(i, res.z());
		}


		// Build matrix A
//		long selected = Arrays.stream(mesh.getVertices()).filter(vertex -> vertex.hasTag(PsObject.IS_SELECTED)).count();
		PnSparseMatrix A = new PnSparseMatrix(mesh.getNumVertices(), mesh.getNumVertices());
		for (int i = 0; i < mesh.getNumVertices(); i++) {
			A.setEntry(i, i, constrained.contains(i) ? 1 : 0);
		}

		PnSparseMatrix At = PnSparseMatrix.transposeNew(A);

		// Step 3: build a vectors
		PdVector ax = rightMultVector(A, vx.v, null);
		PdVector ay = rightMultVector(A, vy.v, null);
		PdVector az = rightMultVector(A, vz.v, null);

		// Step 4: make LHS matrix
		PnSparseMatrix LtL = multMatrices(m.Lt, m.L, null);
		PnSparseMatrix AtA = multMatrices(At, A, null);
		AtA.multScalar(lambda);

		PnSparseMatrix lhs = addNew(LtL, AtA);
		lhs.validate();

		// Step 4: make RHS matrix
		PdVector rhsx = buildRhs(m.Lt, deltax, At, ax, lambda);
		PdVector rhsy = buildRhs(m.Lt, deltay, At, ay, lambda);
		PdVector rhsz = buildRhs(m.Lt, deltaz, At, az, lambda);

		applyAfterSolving(lhs, rhsx, rhsy, rhsz);
	}

	private void applyAfterSolving(PnSparseMatrix lhs, PdVector rhsx, PdVector rhsy, PdVector rhsz) throws Exception {
		// Step 5: Solve to get vector v
		PdVector vxTilde = new PdVector(mesh.getNumVertices());
		PdVector vyTilde = new PdVector(mesh.getNumVertices());
		PdVector vzTilde = new PdVector(mesh.getNumVertices());

		long factor = PnMumpsSolver.factor(lhs, PnMumpsSolver.Type.SYMMETRIC_POSITIVE_DEFINITE);
		PnMumpsSolver.solve(factor, vxTilde, rhsx);
		PnMumpsSolver.solve(factor, vyTilde, rhsy);
		PnMumpsSolver.solve(factor, vzTilde, rhsz);

		// Step 6: Apply new coordinates in v to mesh
		WVector center = new WVector(mesh.getCenterOfGravity());
		for (int i = 0; i < mesh.getNumVertices(); i++) {
			mesh.setVertex(i, new PdVector(vxTilde.getEntry(i), vyTilde.getEntry(i), vzTilde.getEntry(i)));
		}

		// translate mesh by difference in center of gravity to prevent it from moving
		mesh.translate(center.minus(mesh.getCenterOfGravity()).v);
		mesh.update(mesh);
	}

	private PdVector buildRhs(PnSparseMatrix Lt, PdVector delta, PnSparseMatrix At, PdVector a, double lambda) {
		PdVector LtDelta = rightMultVector(Lt, delta, null);
		PdVector Ata = rightMultVector(At, a, null);
		Ata.multScalar(lambda);

		return PdVector.addNew(LtDelta, Ata);
	}

	private WVector[] computeV() {
		WVector vx = new WVector(new PdVector(mesh.getNumVertices()));
		WVector vy = new WVector(new PdVector(mesh.getNumVertices()));
		WVector vz = new WVector(new PdVector(mesh.getNumVertices()));
		for (int i = 0; i < mesh.getNumVertices(); i++) {
			vx.v.setEntry(i, mesh.getVertex(i).getEntry(0));
			vy.v.setEntry(i, mesh.getVertex(i).getEntry(1));
			vz.v.setEntry(i, mesh.getVertex(i).getEntry(2));
		}

		return new WVector[]{vx, vy, vz};
	}


	//

	private static WMatrix computeA(WVector N, WVector p1, WVector p2, WVector p3) {
		WVector e1 = p3.minus(p2);
		WVector e2 = p1.minus(p3);
		WVector e3 = p2.minus(p1);

		double scalar =  1.0 / (2.0 * area(p1, p2, p3));

		WMatrix A = WMatrix.columns(N.cross(e1), N.cross(e2), N.cross(e3));
		A = A.mult(scalar);
		return A;
	}

	public static void main(String[] args) {
		WVector a = new WVector(-4, 0, 0);
		WVector b = new WVector(0, 0, 3);
		WVector c = new WVector(0, 3, 0);

		WVector N = new WVector(-1, 1, 1);
		System.out.println(computeA(N, a, b, c));
	}

	public static double area(WVector a, WVector b, WVector c) {
		WVector ab = b.minus(a);
		WVector ac = c.minus(a);
		WVector cross = ab.cross(ac);

		return 0.5 * cross.length();
	}
}

class Matrices {
	PnSparseMatrix G, Gt, Mv, S, L, Lt, M, Mi;

	public Matrices(PnSparseMatrix g, PnSparseMatrix gt, PnSparseMatrix mv, PnSparseMatrix s, PnSparseMatrix l, PnSparseMatrix lt, PnSparseMatrix m, PnSparseMatrix mi) {
		G = g;
		Gt = gt;
		Mv = mv;
		S = s;
		L = l;
		Lt = lt;
		M = m;
		Mi = mi;
	}
}