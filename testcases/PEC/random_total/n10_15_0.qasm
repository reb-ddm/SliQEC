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
h q[0];
t q[0];
t q[9];
t q[3];
ccx q[8], q[5], q[6];
t q[3];
h q[9];
t q[3];
ccx q[8], q[5], q[2];
s q[5];
h q[9];
h q[9];
s q[5];
t q[4];
t q[5];
cx q[6], q[9];
s q[0];
s q[9];
h q[7];
s q[5];
h q[7];
t q[9];
s q[1];
ccx q[6], q[4], q[2];
ccx q[2], q[1], q[9];
ccx q[1], q[8], q[6];
cx q[4], q[9];
t q[0];
cx q[4], q[9];
ccx q[8], q[9], q[5];
