#pragma once

#ifndef  PSRELIABLITY_CONFIG_H
#define PSRELIABLITY_CONFIG_H

// 元件数目统计
struct Statistic
{
	// 元件总数
	int numComp;
	// 发电机数
	int numGen;
	// 输电线路数
	int numLine;
	// 变压器数
	int numTran;
	// 负荷数
	int numLoad;
};

// 计算配置
struct EvaluateParam
{
	// 计算模式
	int mode;
	// 计算步长
	int step;
	// 计算尺度(负荷的时间尺度)
	int duration;
};

// 总运行配置
struct ReliablityEvaluateConfig
{
	Statistic systemInfo;
	EvaluateParam evaluateInfo;
};

#endif