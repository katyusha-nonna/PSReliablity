#pragma once

#ifndef  PSRELIABLITY_CONFIG_H
#define PSRELIABLITY_CONFIG_H

// Ԫ����Ŀͳ��
struct Statistic
{
	// Ԫ������
	int numComp;
	// �������
	int numGen;
	// �����·��
	int numLine;
	// ��ѹ����
	int numTran;
	// ������
	int numLoad;
};

// ��������
struct EvaluateParam
{
	// ����ģʽ
	int mode;
	// ���㲽��
	int step;
	// ����߶�(���ɵ�ʱ��߶�)
	int duration;
};

// ����������
struct ReliablityEvaluateConfig
{
	Statistic systemInfo;
	EvaluateParam evaluateInfo;
};

#endif