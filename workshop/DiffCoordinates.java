package workshop;

import jv.geom.PgElementSet;
import jv.project.PgGeometry;
import jv.vecmath.PdVector;
import jv.vecmath.PiVector;
import jvx.project.PjWorkshop;

import java.awt.*;
import java.util.Queue;
import java.util.*;

public class DiffCoordinates extends PjWorkshop {

	PgElementSet m_geom;
	PgElementSet m_geomSave;

	public DiffCoordinates() {
		super("Differential Coordinates");
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


	public void drawGradients() {
		// TODO
	}
}