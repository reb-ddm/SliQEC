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
cx q[4], q[5];
cx q[2], q[5];
cx q[1], q[6];
t q[4];
t q[1];
ccx q[2], q[7], q[6];
s q[3];
h q[7];
s q[0];
ccx q[2], q[5], q[6];
h q[4];
s q[2];
h q[8];
t q[7];
ccx q[1], q[5], q[0];
h q[5];
h q[1];
cx q[7], q[6];
ccx q[5], q[4], q[7];
ccx q[0], q[3], q[5];
t q[3];
ccx q[0], q[7], q[5];
h q[8];
h q[0];
h q[2];
t q[7];
cx q[8], q[1];
cx q[1], q[0];
tdg q[0];
sdg q[2];
cx q[2], q[3];
t q[7];
s q[4];
cx q[4], q[7];
cx q[8], q[4];
h q[7];
cx q[9], q[1];
