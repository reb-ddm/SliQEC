OPENQASM 2.0;
include "qelib1.inc";
qreg q[33];
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
h q[20];
h q[21];
h q[22];
h q[23];
h q[24];
h q[25];
h q[26];
h q[27];
h q[28];
h q[29];
h q[27];
h q[25];
s q[19];
t q[3];
cx q[8], q[5];
s q[14];
s q[21];
ccx q[3], q[0], q[13];
h q[23];
h q[11];
ccx q[11], q[7], q[15];
ccx q[4], q[25], q[8];
t q[27];
t q[9];
s q[12];
s q[5];
ccx q[26], q[11], q[19];
h q[2];
t q[11];
ccx q[6], q[14], q[17];
ccx q[9], q[29], q[2];
cx q[21], q[4];
t q[13];
ccx q[21], q[7], q[26];
ccx q[11], q[12], q[21];
h q[22];
ccx q[27], q[0], q[3];
ccx q[26], q[13], q[22];
ccx q[4], q[15], q[8];
cx q[26], q[14];
s q[7];
t q[22];
t q[4];
t q[2];
cx q[1], q[15];
ccx q[16], q[18], q[22];
s q[19];
t q[5];
s q[15];
s q[13];
t q[7];
cx q[3], q[11];
t q[15];
t q[14];
h q[23];
h q[12];
ccx q[28], q[29], q[6];
t q[14];
ccx q[14], q[18], q[27];
t q[12];
h q[6];
t q[13];
h q[20];
cx q[21], q[5];
s q[3];
cx q[8], q[27];
h q[23];
h q[28];
h q[21];
t q[19];
ccx q[8], q[12], q[14];
t q[10];
h q[14];
s q[7];
h q[29];
ccx q[27], q[23], q[5];
h q[1];
s q[18];
s q[19];
h q[18];
h q[19];
t q[7];
cx q[29], q[20];
cx q[17], q[19];
ccx q[28], q[11], q[29];
t q[8];
t q[11];
ccx q[0], q[13], q[21];
t q[23];
t q[0];
ccx q[0], q[19], q[5];
h q[9];
t q[24];
s q[23];
ccx q[22], q[9], q[6];
ccx q[24], q[7], q[27];
h q[27];
s q[19];
s q[13];
h q[5];
tdg q[1];
cx q[0], q[1];
tdg q[0];
tdg q[3];
cx q[2], q[3];
tdg q[2];
tdg q[4];
cx q[5], q[4];
sdg q[4];
cx q[7], q[6];
tdg q[6];
x q[8];
sdg q[9];
cx q[12], q[11];
tdg q[11];
cx q[12], q[11];
sdg q[11];
tdg q[11];
tdg q[14];
cx q[13], q[14];
sdg q[13];
h q[23];
t q[15];
h q[23];
h q[15];
s q[22];
t q[24];
cx q[18], q[28];
ccx q[20], q[21], q[19];
h q[22];
t q[20];
cx q[22], q[20];
h q[29];
t q[27];
cx q[26], q[28];
cx q[20], q[19];
cx q[30], q[12];
cx q[31], q[10];
cx q[32], q[28];
