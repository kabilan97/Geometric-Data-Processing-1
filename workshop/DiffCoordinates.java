package workshop;

import jv.geom.PgElementSet;
import jv.object.PsDebug;
import jv.project.PgGeometry;
import jv.vecmath.PdMatrix;
import jv.vecmath.PiVector;
import jvx.numeric.PnSparseMatrix;
import jvx.project.PjWorkshop;

import java.util.HashMap;

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


	public void drawGradients() {
		mesh.makeElementNormals();
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

		PnSparseMatrix Gt = G.transposeNew();
		PnSparseMatrix S = PnSparseMatrix.multMatrices(Gt, PnSparseMatrix.multMatrices(Mv, G, null), null);
		PnSparseMatrix Mi = new PnSparseMatrix(mesh.getNumVertices(), mesh.getNumVertices());
		for (int i = 0; i < mesh.getNumVertices(); i++) {
			Mi.setEntry(i, i, 1.0 / M.getEntry(i, i));
		}

		PnSparseMatrix L = PnSparseMatrix.multMatrices(Mi, S, null);
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