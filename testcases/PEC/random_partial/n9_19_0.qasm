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
h q[8];
cx q[7], q[6];
s q[0];
h q[4];
t q[8];
h q[5];
t q[5];
cx q[1], q[2];
h q[1];
cx q[3], q[0];
t q[3];
t q[1];
cx q[6], q[4];
h q[0];
s q[4];
ccx q[0], q[1], q[6];
ccx q[4], q[2], q[1];
s q[3];
t q[1];
h q[3];
cx q[1], q[0];
cx q[2], q[5];
cx q[6], q[5];
cx q[2], q[0];
s q[6];
s q[1];
t q[8];
h q[4];
cx q[1], q[0];
tdg q[0];
cx q[2], q[3];
sdg q[3];
cx q[2], q[3];
sdg q[2];
h q[4];
t q[6];
t q[5];
s q[8];
h q[5];
