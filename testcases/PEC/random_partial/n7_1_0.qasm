OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
s q[1];
h q[1];
h q[0];
ccx q[4], q[0], q[6];
ccx q[3], q[6], q[5];
cx q[0], q[5];
ccx q[6], q[5], q[4];
t q[5];
ccx q[1], q[6], q[4];
t q[4];
t q[5];
cx q[4], q[6];
h q[3];
ccx q[4], q[1], q[3];
t q[4];
t q[5];
h q[6];
t q[4];
t q[5];
ccx q[2], q[3], q[5];
h q[4];
tdg q[0];
t q[3];
h q[5];
ccx q[3], q[5], q[6];
ccx q[3], q[6], q[5];
