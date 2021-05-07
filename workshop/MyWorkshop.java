package workshop;

import java.awt.Color;
import java.util.*;


import jv.geom.PgElementSet;
import jv.object.PsDebug;
import jv.project.PgGeometry;
import jv.vecmath.PdVector;
import jv.vecmath.PdMatrix;
import jv.vecmath.PiVector;
import jvx.project.PjWorkshop;
import jvx.numeric.PnStiffDiriConforming;

import util.Util;

public class MyWorkshop extends PjWorkshop {

	PgElementSet m_geom;
	PgElementSet m_geomSave;
	
	public MyWorkshop() {
		super("My Workshop");
		init();
	}
	
	@Override
	public void setGeometry(PgGeometry geom) {
		super.setGeometry(geom);
		m_geom 		= (PgElementSet)super.m_geom;
		m_geomSave 	= (PgElementSet)super.m_geomSave;
	}
	
	public void init() {		
		super.init();
	}
	
	public void makeRandomElementColors() {
		//assure that the color array is allocated
		m_geom.assureElementColors();
		
		Random rand = new Random();
		Color randomColor;
		
		int noe = m_geom.getNumElements();
		for(int i=0; i<noe; i++){
			randomColor = Color.getHSBColor(rand.nextFloat(), 1.0f, 1.0f);//new Color(rand.nextFloat(), rand.nextFloat(), rand.nextFloat());
			m_geom.setElementColor(i, randomColor);
		}
		m_geom.showElementColorFromVertices(false);
		m_geom.showElementColors(true);	
		m_geom.showSmoothElementColors(false);
	}
	
	
	public void makeRandomVertexColors() {
		//assure that the color array is allocated
		m_geom.assureVertexColors();
		
		Random rand = new Random();
		Color randomColor;
		
		int nov = m_geom.getNumVertices();
		for(int i=0; i<nov; i++){
			randomColor = Color.getHSBColor(rand.nextFloat(), 1.0f, 1.0f);
			m_geom.setVertexColor(i, randomColor);
		}
		
		m_geom.showElementColors(true);	
		m_geom.showVertexColors(true);
		m_geom.showElementColorFromVertices(true);	
		m_geom.showSmoothElementColors(true);
		
	}
	
	
	public void setXOff(double xOff) {
		int nov = m_geom.getNumVertices();
		PdVector v = new PdVector(3);
		// the double array is v.m_data 
		for (int i=0; i<nov; i++) {
			v.copyArray(m_geomSave.getVertex(i));
			v.setEntry(0, v.getEntry(0)+xOff);
			m_geom.setVertex(i, v);
		}
	}

	public String modelInfo() {
		return "Vertices: " + m_geom.getVertices().length + ", Triangles: " + m_geom.getElements().length;
	}

	/**
	 * Compute the volume using signed volume of tetrahedrons
	 */
	public double computeVolume() {
		PdVector[] vertices = m_geom.getVertices();

		PdVector center = average(vertices);

		// construct tetrahedron of each triangle using it's three vertices and the center point
		double volume = Arrays.stream(m_geom.getElements()).mapToDouble(indices -> {
			PdVector p0 = vertices[indices.getEntry(0)];
			PdVector p1 = vertices[indices.getEntry(1)];
			PdVector p2 = vertices[indices.getEntry(2)];
			return signedVolume(p0, p1, p2, center);
		}).sum();



		return Math.abs(volume);
	}

	// get average point from a list of vectors
	private PdVector average(PdVector[] v) {
		PdVector sum = Arrays.stream(v).reduce(PdVector::addNew).get();
		sum.multScalar(1.0 / v.length);
		return sum;
	}

	// https://stackoverflow.com/a/9866530
	private double signedVolume(PdVector a, PdVector b, PdVector c, PdVector d) {
		PdVector v1 = PdVector.subNew(a, d);
		PdVector v2 = PdVector.subNew(b, d);
		PdVector v3 = PdVector.subNew(c, d);

		return PdVector.dot(v1, PdVector.crossNew(v2, v3)) / 6.0;
	}


	public int computeGenus() {
		/** assuming the numbers are defined since the tool itself does show them,
		 * the genus will look as following
		 * defined by the euler formula: (E - V - F - b)/2 + 1
		 **/
		return (m_geom.getNumEdges() - m_geom.getNumVertices() - m_geom.getNumElements()- m_geom.getNumBoundaries())/2 +1;
	}


	public int computeConnectedComponents() {
		m_geom.makeNeighbour();

		Set<Integer> discovered = new HashSet<>();//BFS(0);
		int count = 0;
		for (int i = 0; i < m_geom.getNumElements(); i++) {
			if(!discovered.contains(i)){
				discovered.addAll(BFS(i));
				count++;
			}
		}


		return count;
	}

	private Set<Integer> BFS(int root) {
		Set<Integer> discovered = new HashSet<>();
		Queue<Integer> queue = new ArrayDeque<Integer>();
		discovered.add(root);
		queue.add(root);
		while (!queue.isEmpty()) {
			PiVector neighbours = m_geom.getNeighbour(queue.poll());
			if (neighbours != null) {
				Set<Integer> tmp = new HashSet<>();

				int[] entries = neighbours.getEntries();
				for (int entry : entries) {
					tmp.add(entry);
				}
				tmp.removeAll(discovered);
				queue.addAll(tmp);
				discovered.addAll(tmp);
			}
		}
		return discovered;
	}



}