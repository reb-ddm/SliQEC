OPENQASM 2.0;
include "qelib1.inc";
qreg q[9];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
t q[2];
cx q[7], q[3];
s q[7];
ccx q[6], q[5], q[2];
s q[1];
ccx q[7], q[6], q[3];
ccx q[4], q[2], q[0];
h q[2];
s q[7];
ccx q[2], q[0], q[3];
s q[7];
s q[0];
s q[2];
s q[0];
t q[5];
cx q[4], q[6];
t q[5];
t q[1];
t q[6];
h q[7];
ccx q[5], q[4], q[2];
cx q[1], q[5];
h q[0];
s q[5];
tdg q[0];
x q[2];
cx q[7], q[5];
t q[7];
h q[5];
t q[6];
cx q[8], q[3];