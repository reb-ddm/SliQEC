OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
h q[0];
h q[1];
h q[2];
t q[0];
s q[1];
t q[2];
ccx q[0], q[1], q[2];
h q[1];
t q[1];
s q[1];
t q[1];
h q[0];
t q[1];
s q[2];
