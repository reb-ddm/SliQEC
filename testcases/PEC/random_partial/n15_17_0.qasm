OPENQASM 2.0;
include "qelib1.inc";
qreg q[15];
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
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
t q[4];
h q[10];
t q[11];
t q[7];
ccx q[10], q[5], q[8];
s q[3];
cx q[3], q[5];
ccx q[0], q[1], q[13];
h q[14];
cx q[8], q[12];
s q[10];
cx q[1], q[6];
cx q[3], q[0];
ccx q[6], q[2], q[8];
ccx q[5], q[12], q[13];
t q[10];
s q[13];
t q[13];
h q[12];
ccx q[10], q[1], q[5];
cx q[0], q[2];
t q[1];
cx q[12], q[7];
s q[9];
t q[6];
s q[13];
ccx q[14], q[0], q[2];
t q[11];
s q[12];
cx q[0], q[3];
cx q[9], q[1];
h q[9];
s q[9];
s q[12];
ccx q[8], q[14], q[5];
s q[11];
s q[0];
cx q[8], q[12];
ccx q[7], q[8], q[5];
ccx q[14], q[5], q[12];
cx q[1], q[6];
cx q[12], q[7];
ccx q[13], q[7], q[6];
ccx q[2], q[13], q[6];
cx q[0], q[6];
z q[0];
z q[1];
cx q[2], q[3];
tdg q[2];
tdg q[5];
cx q[4], q[5];
sdg q[4];
cx q[8], q[7];
cx q[9], q[12];
cx q[11], q[12];
h q[12];
h q[14];
h q[7];
t q[12];
h q[11];
