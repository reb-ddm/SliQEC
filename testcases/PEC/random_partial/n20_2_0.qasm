OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
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
h q[15];
h q[16];
h q[17];
h q[18];
h q[19];
ccx q[16], q[4], q[13];
ccx q[6], q[14], q[13];
s q[11];
t q[19];
t q[18];
cx q[17], q[7];
cx q[10], q[1];
s q[5];
h q[12];
cx q[7], q[0];
ccx q[13], q[2], q[15];
t q[5];
h q[14];
t q[13];
t q[10];
h q[0];
ccx q[16], q[19], q[14];
cx q[1], q[2];
h q[4];
ccx q[2], q[12], q[8];
t q[4];
t q[15];
cx q[14], q[12];
t q[18];
t q[11];
ccx q[5], q[12], q[2];
ccx q[15], q[10], q[19];
h q[2];
s q[2];
t q[8];
s q[18];
s q[13];
h q[19];
cx q[17], q[0];
s q[10];
s q[1];
h q[12];
ccx q[1], q[11], q[2];
ccx q[19], q[5], q[2];
t q[19];
ccx q[14], q[0], q[16];
ccx q[3], q[5], q[0];
s q[9];
h q[10];
s q[11];
s q[1];
ccx q[12], q[19], q[3];
h q[12];
ccx q[16], q[19], q[3];
h q[16];
ccx q[16], q[14], q[2];
t q[7];
cx q[19], q[8];
ccx q[19], q[1], q[0];
s q[12];
ccx q[6], q[14], q[4];
ccx q[10], q[9], q[18];
h q[0];
h q[7];
ccx q[16], q[3], q[13];
x q[0];
sdg q[2];
tdg q[2];
cx q[1], q[2];
tdg q[1];
cx q[3], q[4];
tdg q[4];
cx q[3], q[4];
sdg q[4];
tdg q[3];
tdg q[5];
cx q[5], q[6];
tdg q[5];
y q[7];
cx q[9], q[8];
tdg q[9];
tdg q[8];
s q[16];
t q[11];
cx q[16], q[12];
ccx q[15], q[18], q[13];
cx q[16], q[13];
s q[15];
ccx q[18], q[17], q[14];
h q[12];
s q[15];
s q[10];
