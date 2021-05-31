package workshop;

import jv.number.PuDouble;
import jv.object.PsDebug;
import jv.object.PsDialog;
import jv.object.PsUpdateIf;
import jvx.project.PjWorkshop_IP;

import java.awt.*;
import java.io.PrintWriter;
import java.io.StringWriter;

public class DiffCoordinates_IP extends PjWorkshop_IP {

	protected Button modifyMesh;
	protected Button m_bMakeRandomVertexColors;
	protected PuDouble m_xOff;

	protected TextArea matrixInput;
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

		add(new Label("Matrix input"));
		matrixInput = new TextArea("1 0 0\n0 1 0\n0 0 1", 3, 12);
		add(matrixInput);

//		add(new Label("Genus: "));
		modifyMesh = new Button("Update mesh");
		modifyMesh.addActionListener(e -> {
			PsDebug.message("Running");
			try {
				String[] matrix = matrixInput.getText().split("\\s");
				WMatrix m = new WMatrix(new double[][]{
						{d(matrix[0]), d(matrix[1]), d(matrix[2])},
						{d(matrix[3]), d(matrix[4]), d(matrix[5])},
						{d(matrix[6]), d(matrix[7]), d(matrix[8])}
				});

				diffCoordinates.deformMesh(m);
			} catch (Exception ex) {
				PsDebug.message("Exception occurred");
				StringWriter sw = new StringWriter();
				PrintWriter pw = new PrintWriter(sw);
				ex.printStackTrace(pw);

				PsDebug.message(sw.toString());
			}
			PsDebug.message("Done");
		});
		add(modifyMesh);

		validate();
	}

	private static double d(String s) {
		return Double.parseDouble(s);
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
