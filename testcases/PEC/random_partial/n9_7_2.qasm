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
s q[7];
s q[3];
t q[0];
h q[7];
t q[8];
cx q[0], q[8];
ccx q[0], q[2], q[4];
h q[3];
cx q[8], q[3];
h q[3];
cx q[0], q[2];
h q[2];
cx q[3], q[7];
h q[5];
ccx q[6], q[7], q[8];
ccx q[3], q[8], q[0];
t q[6];
cx q[6], q[2];
t q[4];
ccx q[4], q[0], q[1];
s q[3];
cx q[0], q[1];
cx q[5], q[1];
cx q[0], q[5];
cx q[7], q[0];
t q[4];
s q[2];
cx q[0], q[1];
sdg q[1];
tdg q[0];
cx q[2], q[3];
sdg q[3];
cx q[2], q[3];
sdg q[2];
cx q[8], q[4];
t q[8];
t q[6];
cx q[7], q[8];
cx q[7], q[8];
cx q[9], q[5];
