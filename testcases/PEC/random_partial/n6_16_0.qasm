OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
s q[2];
h q[5];
h q[5];
ccx q[3], q[0], q[5];
s q[3];
ccx q[3], q[2], q[0];
s q[5];
ccx q[3], q[4], q[0];
h q[5];
h q[0];
s q[3];
cx q[4], q[1];
s q[5];
ccx q[0], q[5], q[4];
cx q[3], q[2];
s q[3];
h q[5];
t q[5];
cx q[1], q[0];
sdg q[1];
tdg q[0];
cx q[1], q[0];
sdg q[0];
ccx q[3], q[5], q[4];
cx q[5], q[3];
s q[4];
