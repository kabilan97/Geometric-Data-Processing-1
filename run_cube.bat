@echo off
set jv_jar=jars/javaview.jar;jars/jvx.jar;jars/vgpapp.jar;jars\Jama-1.0.3.jar;.
start javaw -cp %jv_jar% -Djava.library.path="dll" -Xmx1024m javaview model="models/custom/cube_hq_unmoved.obj;models/custom/cube_hq_4deg_move0.2.obj" codebase=. archive.dev=show %*