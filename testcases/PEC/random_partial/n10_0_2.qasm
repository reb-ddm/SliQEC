OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
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
h q[2];
cx q[6], q[4];
s q[6];
t q[2];
h q[1];
cx q[0], q[4];
s q[2];
cx q[2], q[3];
h q[6];
s q[3];
s q[7];
t q[9];
ccx q[1], q[0], q[6];
h q[6];
t q[4];
t q[8];
ccx q[8], q[1], q[3];
t q[1];
t q[0];
s q[8];
cx q[1], q[9];
ccx q[1], q[4], q[0];
cx q[0], q[8];
cx q[6], q[1];
t q[3];
cx q[2], q[0];
h q[3];
ccx q[2], q[8], q[4];
ccx q[7], q[3], q[2];
cx q[7], q[1];
z q[0];
tdg q[1];
sdg q[3];
cx q[4], q[3];
tdg q[3];
ccx q[5], q[9], q[8];
t q[6];
cx q[9], q[6];
s q[8];
s q[9];
cx q[10], q[0];
cx q[11], q[5];
