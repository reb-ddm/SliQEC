OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];
h q[9];
t q[1];
s q[1];
t q[5];
cx q[2], q[0];
ccx q[7], q[0], q[1];
cx q[3], q[4];
t q[7];
ccx q[2], q[3], q[0];
s q[2];
t q[8];
t q[1];
ccx q[7], q[2], q[0];
cx q[6], q[8];
t q[5];
cx q[4], q[2];
ccx q[0], q[7], q[1];
t q[0];
ccx q[4], q[2], q[3];
cx q[5], q[4];
t q[9];
ccx q[2], q[4], q[6];
cx q[1], q[0];
ccx q[3], q[5], q[2];
s q[3];
cx q[6], q[9];
h q[6];
ccx q[6], q[0], q[2];
cx q[1], q[4];
s q[7];
ccx q[7], q[8], q[0];
