package workshop;

import jv.number.PuDouble;
import jv.object.PsDebug;
import jv.object.PsDialog;
import jv.object.PsUpdateIf;
import jvx.project.PjWorkshop_IP;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.PrintWriter;
import java.io.StringWriter;

public class DiffCoordinates_IP extends PjWorkshop_IP {

	protected Button draw_gradients;
	protected Button m_bMakeRandomVertexColors;
	protected PuDouble m_xOff;

	DiffCoordinates diffCoordinates;

	public DiffCoordinates_IP() {
		super();
		if(getClass() == DiffCoordinates_IP.class)
			init();
	}

	
	public String getNotice() {
		return "Differential coordinates";
	}
	
	public void setParent(PsUpdateIf parent) {
		super.setParent(parent);
		diffCoordinates = (DiffCoordinates)parent;


		addSubTitle("Task 1");
//		add(new Label("Genus: "));
		draw_gradients = new Button("Draw gradients");
		draw_gradients.addActionListener(e -> {
			PsDebug.message("Running");
			try {
				diffCoordinates.drawGradients();
			} catch (Exception ex) {
				PsDebug.message("Exception occurred");
				StringWriter sw = new StringWriter();
				PrintWriter pw = new PrintWriter(sw);
				ex.printStackTrace(pw);

				PsDebug.message(sw.toString());
			}
			PsDebug.message("Done");
		});
		add(draw_gradients);

		validate();
	}

	public void init() {
		super.init();
		setTitle("Differential Coordinates");
	}
	
	public boolean update(Object event) {
		if (event == m_xOff) {

			return true;
		} else
			return super.update(event);
	}
	/**
	 * Get information which bottom buttons a dialog should create
	 * when showing this info panel.
	 */
	protected int getDialogButtons()		{
		return PsDialog.BUTTON_OK;
	}
}
