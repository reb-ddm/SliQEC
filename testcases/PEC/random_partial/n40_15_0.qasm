OPENQASM 2.0;
include "qelib1.inc";
qreg q[40];
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
h q[30];
h q[31];
h q[32];
h q[33];
h q[34];
h q[35];
h q[36];
h q[37];
h q[38];
h q[39];
t q[13];
cx q[11], q[24];
t q[31];
t q[11];
t q[18];
s q[28];
s q[9];
ccx q[11], q[31], q[15];
cx q[6], q[16];
ccx q[0], q[11], q[30];
s q[30];
t q[24];
cx q[19], q[37];
cx q[30], q[31];
h q[34];
t q[37];
ccx q[2], q[4], q[5];
h q[20];
h q[24];
s q[38];
h q[37];
s q[3];
cx q[21], q[38];
h q[32];
ccx q[14], q[13], q[16];
t q[30];
t q[31];
h q[8];
t q[36];
ccx q[11], q[36], q[9];
t q[30];
s q[8];
cx q[24], q[1];
t q[28];
t q[21];
h q[31];
h q[4];
s q[3];
h q[5];
t q[1];
t q[8];
s q[8];
t q[38];
ccx q[35], q[28], q[8];
t q[30];
cx q[13], q[2];
ccx q[16], q[11], q[28];
ccx q[17], q[20], q[3];
t q[6];
s q[36];
cx q[1], q[18];
ccx q[17], q[19], q[39];
h q[29];
ccx q[26], q[19], q[5];
s q[30];
ccx q[4], q[3], q[8];
h q[38];
ccx q[3], q[16], q[31];
ccx q[34], q[26], q[4];
cx q[36], q[27];
cx q[35], q[24];
ccx q[23], q[21], q[2];
cx q[37], q[8];
s q[19];
cx q[4], q[7];
cx q[3], q[17];
t q[38];
s q[2];
ccx q[0], q[31], q[1];
t q[28];
ccx q[28], q[35], q[25];
h q[0];
h q[32];
h q[9];
cx q[34], q[23];
s q[5];
ccx q[17], q[27], q[22];
t q[9];
cx q[0], q[37];
h q[34];
t q[32];
ccx q[33], q[39], q[4];
s q[39];
h q[17];
s q[18];
ccx q[29], q[16], q[14];
ccx q[29], q[8], q[38];
ccx q[38], q[10], q[3];
ccx q[15], q[12], q[17];
t q[22];
cx q[23], q[12];
t q[36];
t q[19];
s q[30];
t q[19];
s q[29];
cx q[4], q[30];
cx q[26], q[10];
cx q[31], q[36];
cx q[38], q[34];
s q[23];
h q[3];
h q[14];
ccx q[25], q[8], q[26];
h q[19];
t q[35];
s q[22];
t q[16];
t q[17];
ccx q[18], q[25], q[11];
ccx q[2], q[8], q[32];
cx q[23], q[33];
t q[32];
s q[4];
t q[2];
cx q[27], q[32];
h q[13];
cx q[25], q[12];
t q[39];
ccx q[28], q[6], q[0];
y q[0];
cx q[1], q[2];
tdg q[2];
sdg q[2];
cx q[1], q[2];
sdg q[1];
x q[3];
cx q[6], q[5];
sdg q[5];
z q[7];
sdg q[9];
cx q[9], q[8];
sdg q[8];
x q[10];
sdg q[11];
cx q[13], q[14];
sdg q[13];
y q[15];
z q[16];
cx q[17], q[18];
z q[19];
s q[26];
ccx q[39], q[31], q[20];
s q[39];
s q[35];
t q[31];
t q[22];
cx q[31], q[29];
ccx q[27], q[22], q[23];
s q[20];
s q[35];
cx q[21], q[27];
s q[30];
t q[36];
ccx q[31], q[26], q[27];
s q[38];
ccx q[22], q[28], q[32];
cx q[39], q[31];
t q[26];
s q[39];
ccx q[31], q[32], q[20];
