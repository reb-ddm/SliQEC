OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
t q[3];
h q[3];
ccx q[4], q[3], q[0];
h q[4];
cx q[2], q[1];
t q[1];
cx q[2], q[0];
h q[0];
t q[0];
t q[2];
ccx q[2], q[0], q[4];
s q[4];
s q[2];
ccx q[0], q[3], q[1];
t q[4];
z q[0];
ccx q[3], q[4], q[2];
cx q[3], q[2];
ccx q[2], q[4], q[3];
