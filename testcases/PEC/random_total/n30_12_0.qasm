OPENQASM 2.0;
include "qelib1.inc";
qreg q[30];
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
h q[10];
ccx q[26], q[4], q[23];
h q[28];
t q[18];
t q[28];
ccx q[27], q[25], q[23];
s q[20];
ccx q[25], q[10], q[28];
ccx q[6], q[26], q[4];
ccx q[14], q[28], q[12];
s q[2];
t q[23];
cx q[27], q[10];
t q[9];
cx q[10], q[28];
s q[20];
h q[25];
cx q[4], q[10];
h q[4];
s q[16];
s q[15];
t q[18];
h q[14];
cx q[17], q[4];
ccx q[27], q[12], q[2];
ccx q[27], q[19], q[7];
cx q[20], q[13];
ccx q[10], q[17], q[7];
h q[26];
h q[8];
ccx q[12], q[22], q[4];
s q[7];
ccx q[13], q[27], q[19];
ccx q[27], q[24], q[5];
s q[21];
h q[23];
s q[29];
h q[18];
h q[7];
ccx q[3], q[2], q[20];
h q[6];
s q[4];
s q[10];
cx q[16], q[29];
cx q[2], q[11];
ccx q[2], q[10], q[20];
cx q[3], q[16];
t q[1];
h q[1];
t q[9];
h q[25];
s q[3];
t q[11];
t q[5];
cx q[21], q[13];
s q[20];
s q[2];
h q[13];
s q[15];
ccx q[16], q[14], q[28];
ccx q[22], q[16], q[29];
s q[24];
h q[22];
h q[0];
cx q[4], q[6];
t q[18];
h q[2];
ccx q[16], q[8], q[25];
s q[8];
t q[13];
h q[23];
cx q[22], q[0];
t q[17];
t q[13];
ccx q[10], q[15], q[22];
t q[1];
h q[15];
s q[3];
h q[6];
t q[7];
h q[20];
ccx q[7], q[29], q[22];
s q[16];
t q[13];
t q[25];
cx q[28], q[23];
cx q[19], q[23];
ccx q[28], q[13], q[6];
s q[21];
h q[22];