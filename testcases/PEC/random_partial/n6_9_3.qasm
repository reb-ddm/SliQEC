OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
t q[4];
t q[1];
h q[5];
cx q[3], q[5];
s q[0];
cx q[2], q[3];
cx q[5], q[3];
cx q[1], q[3];
s q[3];
s q[2];
h q[1];
cx q[2], q[5];
t q[3];
t q[1];
s q[3];
s q[5];
t q[1];
cx q[1], q[0];
cx q[0], q[1];
cx q[4], q[5];
t q[4];
h q[3];
cx q[6], q[4];
