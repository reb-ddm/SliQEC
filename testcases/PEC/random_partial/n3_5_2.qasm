OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
h q[0];
h q[1];
h q[2];
t q[2];
h q[2];
t q[0];
s q[1];
t q[2];
t q[2];
t q[2];
ccx q[0], q[1], q[2];
cx q[2], q[0];
y q[0];
t q[1];
s q[1];
cx q[3], q[0];
