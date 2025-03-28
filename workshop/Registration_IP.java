package workshop;

import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.Label;
import java.awt.List;
import java.awt.Panel;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import jv.geom.PgElementSet;
import jv.number.PuBoolean;
import jv.number.PuInteger;
import jv.object.PsConfig;
import jv.object.PsDebug;
import jv.object.PsDialog;
import jv.object.PsUpdateIf;
import jv.objectGui.PsList;
import jv.project.PgGeometryIf;
import jv.project.PvGeometryIf;
import jv.viewer.PvDisplay;
import jvx.project.PjWorkshop_IP;


/**
 * Info Panel of Workshop for surface registration
 *
 */
public class Registration_IP extends PjWorkshop_IP implements ActionListener{
	private int kDefault = 3;
	private int nDefault = 100;

	protected	List			m_listActive;
	protected	List			m_listPassive;
	protected	Vector			m_geomList;
	protected	Registration	m_registration;
	protected   Button			m_bSetSurfaces;

	protected  PuInteger k;
	protected  PuInteger n;
	protected PuBoolean usePointToPlane;
	protected PuBoolean runInLoop;
	private Label outputLabel;

	/** Constructor */
	public Registration_IP () {
		super();
		if (getClass() == Registration_IP.class)
			init();
	}

	/**
	 * Informational text on the usage of the dialog.
	 * This notice will be displayed if this info panel is shown in a dialog.
	 * The text is split at line breaks into individual lines on the dialog.
	 */
	public String getNotice() {
		return "Select the parts and press the button.";
	}
	
	/** Assign a parent object. */
	public void setParent(PsUpdateIf parent) {
		super.setParent(parent);
		m_registration = (Registration)parent;
		
		Panel pGeometries = new Panel();
		pGeometries.setLayout(new GridLayout(1, 2));

		k = new PuInteger("K value");
		k.setDefValue(kDefault);
		k.setValue(kDefault);
		add(k.getInfoPanel());

		n = new PuInteger("N value");
		n.setDefValue(nDefault);
		n.setValue(nDefault);
		n.setBounds(1, 10000);
		add(n.getInfoPanel());

		usePointToPlane = new PuBoolean("Use Point-to-Plane method");
		add(usePointToPlane.getInfoPanel());

		runInLoop = new PuBoolean("Run until translation < 0.01");
		add(runInLoop.getInfoPanel());

		Panel Passive = new Panel();
		Passive.setLayout(new BorderLayout());
		Panel Active = new Panel();
		Active.setLayout(new BorderLayout());
		Label ActiveLabel = new Label("Surface P");
		Active.add(ActiveLabel, BorderLayout.NORTH);
		m_listActive = new PsList(3, true);
		Active.add(m_listActive, BorderLayout.CENTER);
		pGeometries.add(Active);
		Label PassiveLabel = new Label("Surface Q");
		Passive.add(PassiveLabel, BorderLayout.NORTH);
		m_listPassive = new PsList(3, true);
		Passive.add(m_listPassive, BorderLayout.CENTER);
		pGeometries.add(Passive);
		add(pGeometries);
		
		Panel pSetSurfaces = new Panel(new BorderLayout());
		m_bSetSurfaces = new Button("Set selected surfaces");
		m_bSetSurfaces.addActionListener(this);
		pSetSurfaces.add(m_bSetSurfaces, BorderLayout.CENTER);
		add(pSetSurfaces);

		outputLabel = new Label("None yet");
		add(outputLabel);
		
		updateGeomList();
		validate();
	}
		
	/** Initialisation */
	public void init() {
		super.init();
		setTitle("Surface Registration");
		
	}

	/** Set the list of geometries in the lists to the current state of the display. */
	public void updateGeomList() {
		Vector displays = m_registration.getGeometry().getDisplayList();
		int numDisplays = displays.size();
		m_geomList = new Vector();
		for (int i=0; i<numDisplays; i++) {
			PvDisplay disp =((PvDisplay)displays.elementAt(i));
			PgGeometryIf[] geomList = disp.getGeometries();
			int numGeom = geomList.length;
			for (int j=0; j<numGeom; j++) {
				if (!m_geomList.contains(geomList[j])) {
					//Take just PgElementSets from the list.
					if (geomList[j].getType() == PvGeometryIf.GEOM_ELEMENT_SET)
						m_geomList.addElement(geomList[j]);
				}
			}
		}
		int nog = m_geomList.size();
		m_listActive.removeAll();
		m_listPassive.removeAll();
		for (int i=0; i<nog; i++) {
			String name = ((PgGeometryIf)m_geomList.elementAt(i)).getName();
			m_listPassive.add(name);
			m_listActive.add(name);
		}
	}
	/**
	 * Handle action events fired by buttons etc.
	 */
	public void actionPerformed(ActionEvent event) {
		Object source = event.getSource();
		if (source == m_bSetSurfaces) {
			m_registration.setGeometries((PgElementSet)m_geomList.elementAt(m_listActive.getSelectedIndex()),
			(PgElementSet)m_geomList.elementAt(m_listPassive.getSelectedIndex()));

			if (runInLoop.getState()) {
				for (int i = 0; i < 100; i++) {
					outputLabel.setText("Computing");
					double total = m_registration.run(n.getValue(), k.getValue(), usePointToPlane.getState());
					outputLabel.setText("Done (Iterations: " +  (i+1)  + ", translation = " + total + ") ");

					if (total < 0.01) {
						outputLabel.setText("Done - ALIGNED (Iterations: " +  (i+1)  + ", translation = " + total + ")");
						break;
					}
				}
			} else {
				outputLabel.setText("Computing");
				double total = m_registration.run(n.getValue(), k.getValue(), usePointToPlane.getState());
				outputLabel.setText("Done (translation = " + total + ")");
			}
		}
	}
	/**
	 * Get information which bottom buttons a dialog should create
	 * when showing this info panel.
	 */
	protected int getDialogButtons()		{
		return PsDialog.BUTTON_OK;
	}
}
