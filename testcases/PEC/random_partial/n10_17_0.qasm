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
cx q[2], q[1];
ccx q[4], q[9], q[7];
s q[1];
ccx q[1], q[6], q[9];
t q[2];
ccx q[3], q[6], q[9];
s q[9];
h q[7];
ccx q[3], q[9], q[1];
ccx q[6], q[9], q[1];
t q[8];
h q[6];
h q[2];
s q[1];
h q[6];
cx q[6], q[8];
t q[9];
s q[3];
h q[9];
cx q[9], q[7];
ccx q[5], q[6], q[1];
ccx q[9], q[4], q[5];
h q[2];
h q[0];
t q[9];
ccx q[5], q[1], q[3];
h q[7];
cx q[7], q[3];
h q[8];
s q[4];
tdg q[0];
cx q[1], q[0];
sdg q[0];
cx q[2], q[3];
sdg q[2];
ccx q[5], q[8], q[6];
h q[6];
h q[6];
ccx q[5], q[8], q[7];
t q[8];
