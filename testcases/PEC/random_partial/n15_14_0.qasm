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
s q[5];
t q[6];
ccx q[12], q[0], q[11];
s q[11];
cx q[10], q[11];
cx q[5], q[7];
t q[1];
ccx q[1], q[13], q[10];
t q[11];
s q[8];
ccx q[10], q[11], q[0];
t q[6];
t q[2];
s q[11];
h q[9];
ccx q[7], q[12], q[10];
s q[11];
ccx q[6], q[3], q[8];
ccx q[4], q[2], q[0];
ccx q[13], q[8], q[7];
cx q[0], q[9];
cx q[2], q[8];
cx q[1], q[0];
h q[13];
h q[10];
s q[11];
s q[9];
ccx q[11], q[5], q[4];
cx q[11], q[10];
s q[10];
t q[14];
t q[14];
t q[3];
t q[13];
t q[1];
h q[6];
ccx q[8], q[7], q[13];
ccx q[3], q[9], q[1];
cx q[1], q[10];
t q[7];
ccx q[3], q[12], q[9];
cx q[8], q[5];
cx q[9], q[10];
ccx q[4], q[0], q[1];
t q[2];
sdg q[0];
y q[2];
sdg q[5];
cx q[4], q[5];
sdg q[4];
y q[6];
cx q[8], q[7];
h q[14];
s q[14];
t q[8];
cx q[9], q[13];
ccx q[11], q[14], q[8];
ccx q[10], q[13], q[12];
h q[7];
