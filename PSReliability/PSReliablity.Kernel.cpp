#include "PSReliablity.Kernel.h"

OutageTable* CreateTableForComponent(BaseModel* cp, int start, int length, int duration)
{
	// 存放结果
	auto result = new OutageTable(start);
	result->exactProbability.resize(length);
	result->incrementalFrequency.resize(length);
	result->accumulatedProbability.resize(length);
	result->accumulatedFrequency.resize(length);
	// 填充结果
	// 首先求解强迫停运率q和运行概率p
	if (auto gen = dynamic_cast<GeneratorModel*>(cp))
	{
		// 对于发电机元件
		// 首先计算强迫停运率和运行概率
		double q = gen->lambda / (gen->lambda + gen->miu);
		double p = 1 - q;
		// 然后计算增量频率
		double f = q * gen->miu;
		// 然后生成停运表
		result->exactProbability[0] = p;
		result->exactProbability[length - 1] = q;
		for (int i = 1; i < length; i++)
		{
			result->accumulatedProbability[i] = q;
		}
		result->accumulatedProbability[0] = 1;
		result->incrementalFrequency[0] = -f;
		result->incrementalFrequency[length - 1] = f;
		for (int i = 1; i < length; i++)
		{
			result->accumulatedFrequency[i] = f;
		}
		result->accumulatedFrequency[0] = 0;
	}
	else if (auto load = dynamic_cast<LoadModel*>(cp))
	{
		// 对于负荷元件
		// 统计各负荷值出现次数，生成确切概率和累计概率
		for(int i = 0; i < length; i++)
		{
			double num = std::count(load->loadCurve.begin(), load->loadCurve.end(), i);
			result->exactProbability[i] = num / duration;
		}
		result->accumulatedProbability[length-1] = result->exactProbability[length - 1];
		for(int i = length - 2; i > -1; i--)
		{
			result->accumulatedProbability[i] = result->accumulatedProbability[i + 1] + result->exactProbability[i];
		}
		// 扫描负荷曲线，生成增量频率和累计频率
		auto update = [&](int pre, int cur) {
			if (pre < cur)
			{
				// 说明负荷增加
				for (int i = pre; i < cur; i++)
				{
					result->incrementalFrequency[i] += -1;
				}
			}
			else if (pre > cur)
			{
				// 说明负荷减少
				for (int i = pre; i > cur; i--)
				{
					result->incrementalFrequency[i] += 1;
				}
			}
			else
			{
				// 负荷未发生变化，此时不进行操作
			}
		};
		// 先处理第一个时刻的负荷
		update(*(load->loadCurve.rbegin()), *(load->loadCurve.begin()));
		// 然后处理剩下的
		for (int i = 1; i < load->loadCurve.size(); i++)
		{
			update(load->loadCurve[i - 1], load->loadCurve[i]);
		}
		// 生成累积频率
		result->accumulatedFrequency[length - 1] = result->incrementalFrequency[length - 1];
		for (int i = length - 2; i > -1; i--)
		{
			result->accumulatedFrequency[i] = result->accumulatedFrequency[i + 1] + result->incrementalFrequency[i];
		}
	}
	else if (auto line = dynamic_cast<LineModel*>(cp))
	{
		// 对于线路
		// 首先计算强迫停运率和运行概率
		double q = line->lambda / (line->lambda + line->miu);
		double p = 1 - q;
		// 然后计算增量频率
		double f = q * line->miu;
		// 然后生成停运表
		result->exactProbability[0] = p;
		result->exactProbability[length - 1] = q;
		for (int i = 1; i < length; i++)
		{
			result->accumulatedProbability[i] = q;
		}
		result->accumulatedProbability[0] = 1;
		result->incrementalFrequency[0] = -f;
		result->incrementalFrequency[length - 1] = f;
		for (int i = 1; i < length; i++)
		{
			result->accumulatedFrequency[i] = f;
		}
		result->accumulatedFrequency[0] = 0;
	}
	else if (auto tran = dynamic_cast<TransformerModel*>(cp))
	{
		// 对于变压器
				// 首先计算强迫停运率和运行概率
				double q = tran->lambda / (tran->lambda + tran->miu);
				double p = 1 - q;
				// 然后计算增量频率
				double f = q * tran->miu;
				// 然后生成停运表
				result->exactProbability[0] = p;
				result->exactProbability[length - 1] = q;
				for (int i = 1; i < length; i++)
				{
					result->accumulatedProbability[i] = q;
				}
				result->accumulatedProbability[0] = 1;
				result->incrementalFrequency[0] = -f;
				result->incrementalFrequency[length - 1] = f;
				for (int i = 1; i < length; i++)
				{
					result->accumulatedFrequency[i] = f;
				}
				result->accumulatedFrequency[0] = 0;
	}
	return result;
}

OutageTable ParallelConvolution(OutageTable cp1, OutageTable cp2)
{
	// 最终结果的维度、起始时间
	int length = 0, start = 0;
	length = cp1.accumulatedProbability.size() + cp2.accumulatedProbability.size() -1;
	start = cp1.start + cp2.start;
	// 存放结果
	OutageTable result(start);
	std::vector<double> exactProbability(length);
	std::vector<double> accumulatedProbability(length);
	std::vector<double> incrementalFrequency(length);
	std::vector<double> accumulatedFrequency(length);
	// lambda对象，用途：进行卷积运算
	auto convolution = [&](const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& r) {
		int xSize = x.size(), ySize = y.size();
		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < xSize && j <= i; j++)
			{
				if (i - j < ySize) r[i] += x[j] * y[i - j];
			}
		}			
	};
	// lambda对象，用途：改进卷积运算(针对确切概率参与的卷积，且位于y的位置)
	auto convolution2 = [&](const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& r) {
		int xSize = x.size(), ySize = y.size();
		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < xSize; j++)
			{
				if (i - j < ySize)
				{
					if (i - j < 0) r[i] += x[j];
					else r[i] += x[j] * y[i - j];
				}
			}
		}
	};
	// 计算确切概率 = cp1的确切概率和cp2的确切概率的卷积
	convolution(cp1.exactProbability, cp2.exactProbability, exactProbability);
	// 计算累积概率 = cp1的确切概率和cp2的累积概率的卷积
	convolution2(cp1.exactProbability, cp2.accumulatedProbability, accumulatedProbability);
	// 计算增量频率 = cp1的确切概率和cp2的增量频率的卷积 + cp1的增量频率和cp2的确切概率的卷积
	convolution(cp1.exactProbability, cp2.incrementalFrequency, incrementalFrequency);
	convolution(cp1.incrementalFrequency, cp2.exactProbability, incrementalFrequency);
	// 计算累计频率 = cp1的确切概率和cp2的累计频率的卷积 + cp1的增量频率和cp2的累计概率的卷积
	convolution(cp1.exactProbability, cp2.accumulatedFrequency, accumulatedFrequency);
	convolution2(cp1.incrementalFrequency, cp2.accumulatedProbability, accumulatedFrequency);
	// 返回结果
	result.exactProbability = std::move(exactProbability);
	result.accumulatedProbability = std::move(accumulatedProbability);
	result.incrementalFrequency = std::move(incrementalFrequency);
	result.accumulatedFrequency = std::move(accumulatedFrequency);
	return std::move(result);
}

OutageTable SeriesEquivalence(OutageTable cp1, OutageTable cp2)
{
	// 最终结果的维度、起始时间
	int length = 0, start = 0;
	// 首先判断cp1和cp2是否维度相同(容量相同)
	if (cp1.start == cp2.start)
	{
		// 容量相同，直接计算
		length = cp1.exactProbability.size();
	}
	else
	{
		// 容量不同，不能直接计算，要进行改造
		if (cp1.start < cp2.start)
		{
			length = cp1.exactProbability.size();
			int index = cp2.start - cp1.start;
			std::vector<double> temp;
			temp.assign(cp2.accumulatedProbability.begin() + index, cp2.accumulatedProbability.end());
			cp2.accumulatedProbability = std::move(temp);
			cp2.accumulatedProbability[0] = 1;
			temp.assign(cp2.accumulatedFrequency.begin() + index, cp2.accumulatedFrequency.end());
			cp2.accumulatedFrequency = std::move(temp);
			cp2.accumulatedFrequency[0] = 0;
		}
		else
		{
			length = cp2.exactProbability.size();
			int index = cp1.start - cp2.start;
			std::vector<double> temp;
			temp.assign(cp1.accumulatedProbability.begin() + index, cp1.accumulatedProbability.end());
			cp1.accumulatedProbability = std::move(temp);
			cp1.accumulatedProbability[0] = 1;
			temp.assign(cp1.accumulatedFrequency.begin() + index, cp1.accumulatedFrequency.end());
			cp1.accumulatedFrequency = std::move(temp);
			cp1.accumulatedFrequency[0] = 0;
		}
	}
	start = std::min(cp1.start, cp2.start);
	// 存放结果
	OutageTable result(start);
	std::vector<double> exactProbability(length);
	std::vector<double> accumulatedProbability(length);
	std::vector<double> incrementalFrequency(length);
	std::vector<double> accumulatedFrequency(length);
	// 计算累计概率
	for (int i = 0; i < length; i++)
	{
		accumulatedProbability[i] += cp1.accumulatedProbability[i] + cp2.accumulatedProbability[i] - cp1.accumulatedProbability[i] * cp2.accumulatedProbability[i];
	}
	// 计算累计频率
	for (int i = 0; i < length; i++)
	{
		accumulatedFrequency[i] += cp1.accumulatedFrequency[i] * (1 - cp2.accumulatedProbability[i]);
		accumulatedFrequency[i] += cp2.accumulatedFrequency[i] * (1 - cp1.accumulatedProbability[i]);
	}
	// 计算确切概率和增量频率
	exactProbability[length - 1] = accumulatedProbability[length - 1];
	for (int i = length - 2; i > -1; i--)
	{
		exactProbability[i] = accumulatedProbability[i] - accumulatedProbability[i + 1];
	}
	incrementalFrequency[length - 1] = accumulatedFrequency[length - 1];
	for (int i = length - 2; i > -1; i--)
	{
		incrementalFrequency[i] = accumulatedFrequency[i] - accumulatedFrequency[i + 1];
	}
	// 返回结果
	result.exactProbability = std::move(exactProbability);
	result.accumulatedProbability = std::move(accumulatedProbability);
	result.incrementalFrequency = std::move(incrementalFrequency);
	result.accumulatedFrequency = std::move(accumulatedFrequency);
	return std::move(result);
}

void GenerationSystemEvaluate(std::ifstream& input, std::ofstream& output)
{
	// 声明用于存储元件信息和计算中间结果的变量
	// 表达式树，用于存放计算路径
	ExprTree<OutageTable> systemGraph;
	// 字典，用于存放元件列表
	std::map<std::string, BaseModel* >  components;
	// 运行配置
	ReliablityEvaluateConfig config;

	// 必要的初始化
	{
		auto comp = [](std::string a, std::string b)->bool {
			if (a == "+")
			{
				if (b == "*")
				{
					return false;
				}
				else
				{
					return true;
				}
			}
			else if (a == "*")
			{
				return true;
			}
			else
			{
				return false;
			}
		};
		systemGraph.AddComp(comp);
		systemGraph.AddOperator("+", ParallelConvolution);
		systemGraph.AddOperator("*", SeriesEquivalence);
	}

	// 读取模式和步长以及(负荷的)时间尺度
	input >> config.evaluateInfo.mode >> config.evaluateInfo.step >> config.evaluateInfo.duration;
	// 读取系统统计信息
	input >> config.systemInfo.numComp;
	input >> config.systemInfo.numGen;
	input >> config.systemInfo.numLoad;
	input >> config.systemInfo.numLine;
	input >> config.systemInfo.numTran;
	// 读取发电机
	for (int i = 0; i < config.systemInfo.numGen; i++)
	{
		auto gen = new GeneratorModel();
		input >> gen->name >> gen->ID >> gen->capacity >> gen->lambda >> gen->miu;
		components.try_emplace(gen->ID, gen);
	}
	// 读取负荷
	for (int i = 0; i < config.systemInfo.numLoad; i++)
	{
		auto load = new LoadModel();
		input >> load->name >> load->ID >> load->capacity;
		load->loadCurve.resize(config.evaluateInfo.duration);
		for (int i = 0; i < config.evaluateInfo.duration; i++)
		{
			input >> load->loadCurve[i];
			load->loadCurve[i] /= config.evaluateInfo.step;
		}
		components.try_emplace(load->ID, load);
	}
	// 读取线路
	for (int i = 0; i < config.systemInfo.numLine; i++)
	{
		auto line = new LineModel();
		input >> line->name >> line->ID >> line->capacity >> line->lambda >> line->miu;
		components.try_emplace(line->ID, line);
	}
	// 读取变压器
	for (int i = 0; i < config.systemInfo.numTran; i++)
	{
		auto tran = new TransformerModel();
		input >> tran->name >> tran->ID >> tran->capacity >> tran->lambda >> tran->miu;
		components.try_emplace(tran->ID, tran);
	}
	// 为所有元件形成裕度表并填入
	{
		for (auto& cp : components)
		{
			if (auto gen = dynamic_cast<GeneratorModel*>(cp.second))
			{
				int start = gen->capacity / config.evaluateInfo.step;
				int length = start + 1;
				systemGraph.AddData(cp.first, CreateTableForComponent(cp.second, start, length, config.evaluateInfo.duration));
			}
			else if (auto load = dynamic_cast<LoadModel*>(cp.second))
			{
				int start = 0;
				int length = load->capacity / config.evaluateInfo.step + 1;
				systemGraph.AddData(cp.first, CreateTableForComponent(cp.second, start, length, config.evaluateInfo.duration));
			}
			else if (auto line = dynamic_cast<LineModel* >(cp.second))
			{
				int start = line->capacity / config.evaluateInfo.step;
				int length = start + 1;
				systemGraph.AddData(cp.first, CreateTableForComponent(cp.second, start, length, config.evaluateInfo.duration));
			}
			else if (auto tran = dynamic_cast<TransformerModel*>(cp.second))
			{
				int start = tran->capacity / config.evaluateInfo.step;
				int length = start + 1;
				systemGraph.AddData(cp.first, CreateTableForComponent(cp.second, start, length, config.evaluateInfo.duration));
			}
		}
	}
	// 销毁变量
	for (auto& cp : components)
	{
		delete cp.second;
		cp.second = nullptr;
	}
	// 读取网络拓扑并形成表达式树
	{
		std::string s;
		std::stringstream ss;
		while (s == "")
		{
			getline(input, s);
		}
		ss << s;
		systemGraph.BuildFromInfixExpr(ss);
	}
	auto result = systemGraph.Calcute();
	// 输出结果
	{
		// 输出发电系统停运表/裕度表
		int TableLength = 25;
		int count = 0;
		int length = result.exactProbability.size();
		output << "序号" << std::setw(TableLength) << "裕度" << std::setw(TableLength) << "P*" << std::setw(TableLength) << "F*" << std::setw(TableLength)
			<< "P" << std::endl;
		for (int count = 0; count < length; count++)
		{
			output << count << std::setw(TableLength) << (result.start-count) * config.evaluateInfo.step << std::setw(TableLength) << result.accumulatedProbability[count] << std::setw(TableLength) << result.accumulatedFrequency[count] << std::setw(TableLength)
				<< result.exactProbability[count] << std::endl;
		}
		// 输出可靠性指标
		double lolp = 0, eens = 0, lolf = 0;
		int k = result.start + 1;
		lolp = k < result.accumulatedProbability.size() ? result.accumulatedProbability[k] : 0;
		for (int i = k; i < length; i++)
		{
			eens += result.accumulatedProbability[i];
		}
		eens *= config.evaluateInfo.step*config.evaluateInfo.duration;
		lolf = k < result.accumulatedFrequency.size() ? result.accumulatedFrequency[k] : 0;
		output << std::endl;
		output << "LOLP" << std::setw(TableLength) << lolp << std::endl;
		output << "EENS" << std::setw(TableLength) << eens << std::endl;
		output << "LOLF" << std::setw(TableLength) << lolf << std::endl;
	}
}