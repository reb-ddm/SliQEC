OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
h q[0];
h q[1];
h q[2];
h q[3];
ccx q[1], q[0], q[2];
ccx q[3], q[2], q[0];
cx q[3], q[1];
h q[1];
cx q[1], q[3];
ccx q[0], q[1], q[3];
ccx q[3], q[1], q[0];
ccx q[3], q[0], q[1];
ccx q[0], q[1], q[2];
ccx q[3], q[1], q[2];
s q[1];
ccx q[2], q[0], q[1];
cx q[0], q[1];
tdg q[1];
cx q[0], q[1];
sdg q[0];
h q[2];
cx q[3], q[2];
cx q[4], q[2];
