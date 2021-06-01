package workshop;

import jv.number.PuDouble;
import jv.object.PsDebug;
import jv.object.PsDialog;
import jv.object.PsObject;
import jv.object.PsUpdateIf;
import jvx.project.PjWorkshop_IP;

import java.awt.*;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

public class DiffCoordinates_IP extends PjWorkshop_IP {

	protected Button modifyMesh, modifyMeshLaplacian, markSelected, markConstrained;
	protected Button m_bMakeRandomVertexColors;
	protected PuDouble m_xOff;

	protected TextArea matrixInput;
	DiffCoordinates diffCoordinates;

	Set<Integer> selected = new HashSet<>();
	Set<Integer> constrained = new HashSet<>();

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

		markSelected = new Button("Mark selected");
		markConstrained = new Button("Mark constrained");

		markSelected.addActionListener(e -> {
			selected.clear();
			for (int i = 0; i < diffCoordinates.mesh.getNumVertices(); i++) {
				if (diffCoordinates.mesh.getVertex(i).hasTag(PsObject.IS_SELECTED)) {
					selected.add(i);
				}
			}
			markSelected.setLabel("Mark selected (" + selected.size() + ")");
		});
		markConstrained.addActionListener(e -> {
			constrained.clear();
			for (int i = 0; i < diffCoordinates.mesh.getNumVertices(); i++) {
				if (diffCoordinates.mesh.getVertex(i).hasTag(PsObject.IS_SELECTED)) {
					constrained.add(i);
				}
			}
			markConstrained.setLabel("Mark constrained (" + constrained.size() + ")");
		});

		add(markSelected);
		add(markConstrained);

		modifyMesh = new Button("Update mesh");
		modifyMesh.addActionListener(e -> handleButton(false));

		modifyMeshLaplacian = new Button("Update mesh (Laplacian)");
		modifyMeshLaplacian.addActionListener(e -> handleButton(true));
		add(modifyMesh);
		add(modifyMeshLaplacian);

		validate();
	}

	private void handleButton(boolean isLaplacian) {
		PsDebug.message("Running");
		try {
			String[] matrix = matrixInput.getText().split("\\s");
			WMatrix m = new WMatrix(new double[][]{
					{d(matrix[0]), d(matrix[1]), d(matrix[2])},
					{d(matrix[3]), d(matrix[4]), d(matrix[5])},
					{d(matrix[6]), d(matrix[7]), d(matrix[8])}
			});

			if (isLaplacian) {
				diffCoordinates.deformMesh(m, selected);
			} else {
				diffCoordinates.deformMeshLaplacian(m, selected, constrained);
			}
		} catch (Exception ex) {
			PsDebug.message("Exception occurred");
			StringWriter sw = new StringWriter();
			PrintWriter pw = new PrintWriter(sw);
			ex.printStackTrace(pw);

			PsDebug.message(sw.toString());
		}
		selected.clear();
		constrained.clear();
		markConstrained.setLabel("Mark constrained");
		markSelected.setLabel("Mark selected");
		PsDebug.message("Done");
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
