@echo off
javac -classpath jars\javaview.jar;jars\jvx.jar;jars\Jama-1.0.3.jar;. workshop\*.java
javac -classpath jars\javaview.jar;jars\jvx.jar;. menu\*.java
call run_cube.bat