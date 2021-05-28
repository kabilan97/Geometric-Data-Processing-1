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

import java.util.HashMap;

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
			// Compute matrix G
			PiVector vertexIndices = elements[triangleIndex];

			WVector p1 = new WVector(mesh.getVertex(vertexIndices.getEntry(0)));
			WVector p2 = new WVector(mesh.getVertex(vertexIndices.getEntry(1)));
			WVector p3 = new WVector(mesh.getVertex(vertexIndices.getEntry(2)));

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


//		PsDebug.message(Mv.toString());

		PnSparseMatrix Gt = G.transposeNew();
		PnSparseMatrix S = multMatrices(Gt, multMatrices(Mv, G, null), null);
		PnSparseMatrix Mi = new PnSparseMatrix(mesh.getNumVertices(), mesh.getNumVertices());
		for (int i = 0; i < mesh.getNumVertices(); i++) {
			Mi.setEntry(i, i, 1.0 / M.getEntry(i, i));
		}

		PnSparseMatrix L = multMatrices(Mi, S, null);

		return new Matrices(G, Gt, Mv, S, L, M, Mi);
	}

	public void drawGradients() throws Exception {
		mesh.makeElementNormals();

		Matrices m = computeMatrices();

		WMatrix userSuppliedMatrix = new WMatrix(new double[][]{
				{1, 0, 0},
				{0, 1, 0},
				{0, 0, 1}
		});

		// Step 1: build vector g
		WVector vx = new WVector(new PdVector(mesh.getNumVertices()));
		WVector vy = new WVector(new PdVector(mesh.getNumVertices()));
		WVector vz = new WVector(new PdVector(mesh.getNumVertices()));
		for (int i = 0; i < mesh.getNumVertices(); i++) {
			vx.v.setEntry(i, mesh.getVertex(i).getEntry(0));
			vy.v.setEntry(i, mesh.getVertex(i).getEntry(1));
			vz.v.setEntry(i, mesh.getVertex(i).getEntry(2));
		}

		WVector gx = new WVector(rightMultVector(m.G, vx.v, null));
		WVector gy = new WVector(rightMultVector(m.G, vy.v, null));
		WVector gz = new WVector(rightMultVector(m.G, vz.v, null));

		// Step 2: build gTilde vectors
		WVector gxTilde = new WVector(new PdVector(gx.v.getSize()));
		WVector gyTilde = new WVector(new PdVector(gy.v.getSize()));
		WVector gzTilde = new WVector(new PdVector(gz.v.getSize()));

		for (int i = 0; i < gx.v.getSize(); i++) {
			WVector a = new WVector(new PdVector(gx.v.getEntry(i), gy.v.getEntry(i), gz.v.getEntry(i)));

			WVector res = userSuppliedMatrix.mult(a).toVector();

			gxTilde.v.setEntry(i, res.x());
			gyTilde.v.setEntry(i, res.y());
			gzTilde.v.setEntry(i, res.z());
		}

		// Step 3: make LHS matrix
		PnSparseMatrix epsilonM = multScalar(m.M, 0.01);
		PnSparseMatrix lhs = addNew(multMatrices(m.Gt, multMatrices(m.Mv, m.G, null), null), epsilonM);
		lhs.validate();

		// Step 4: make RHS matrix
		PdVector rhsx = rightMultVector(m.Gt, rightMultVector(m.Mv, gxTilde.v, null), null);
		PdVector rhsy = rightMultVector(m.Gt, rightMultVector(m.Mv, gyTilde.v, null), null);
		PdVector rhsz = rightMultVector(m.Gt, rightMultVector(m.Mv, gzTilde.v, null), null);

		// Step 5: Solve to get vector v
		PdVector vxTilde = new PdVector(mesh.getNumVertices());
		PdVector vyTilde = new PdVector(mesh.getNumVertices());
		PdVector vzTilde = new PdVector(mesh.getNumVertices());


//		long factor = PnMumpsSolver.factor(lhs, PnMumpsSolver.Type.GENERAL_SYMMETRIC);
		PnMumpsSolver.solve(lhs, vxTilde, rhsx, PnMumpsSolver.Type.SYMMETRIC_POSITIVE_DEFINITE);
		PnMumpsSolver.solve(lhs, vyTilde, rhsy, PnMumpsSolver.Type.SYMMETRIC_POSITIVE_DEFINITE);
		PnMumpsSolver.solve(lhs, vzTilde, rhsz, PnMumpsSolver.Type.SYMMETRIC_POSITIVE_DEFINITE);

		// Step 6: Apply new coordinates in v to mesh
//		PdVector lhsVx = rightMultVector(lhs, vxTilde, null);
//		PdVector lhsVy = rightMultVector(lhs, vyTilde, null);
//		PdVector lhsVz = rightMultVector(lhs, vzTilde, null);

//		PsDebug.message("X");
//		PsDebug.message(new WVector(lhsVx).toString());
//		PsDebug.message(new WVector(rhsx).toString());
//
//		PsDebug.message("Y");
//		PsDebug.message(new WVector(lhsVy).toString());
//		PsDebug.message(new WVector(rhsy).toString());
//
//		PsDebug.message("Z");
//		PsDebug.message(new WVector(lhsVz).toString());
//		PsDebug.message(new WVector(rhsz).toString());
//
//		PsDebug.message(new WVector(vx.v).toString());
//
//		PsDebug.message(new WVector(vxTilde).toString());
//		PsDebug.message(new WVector(vyTilde).toString());
//		PsDebug.message(new WVector(vzTilde).toString());

		for (int i = 0; i < mesh.getNumVertices(); i++) {
			if (!mesh.getVertex(i).hasTag(PsObject.IS_SELECTED)) {
				continue;
			}

			mesh.setVertex(i, new PdVector(vxTilde.getEntry(i), vyTilde.getEntry(i), vzTilde.getEntry(i)));
		}
		mesh.update(mesh);

	}


	//

	private static WMatrix computeA(WVector N, WVector p1, WVector p2, WVector p3) {
		WVector e1 = p3.minus(p2);
		WVector e2 = p1.minus(p3);
		WVector e3 = p2.minus(p1);

		double scalar =  (1.0 / 2.0) * area(p1, p2, p3);

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
	PnSparseMatrix G, Gt, Mv, S, L, M, Mi;

	public Matrices(PnSparseMatrix g, PnSparseMatrix gt, PnSparseMatrix mv, PnSparseMatrix s, PnSparseMatrix l, PnSparseMatrix m, PnSparseMatrix mi) {
		G = g;
		Gt = gt;
		Mv = mv;
		S = s;
		L = l;
		M = m;
		Mi = mi;
	}
}