# comp_rev函数用于生成互补翻转序列
def comp_rev(seq):
    """
      1. 先将序列反转；
      2. 遍历反转后的序列，将每个碱基替换为其互补碱基：
         A -> T, T -> A, C -> G, G -> C.
    """
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # 使用reversed函数对序列进行反转，并生成互补序列
    return "".join(comp[c] for c in reversed(seq)) # 其中的"".的意思是用空字符串使之相连

def main():
    reference = input("Please input the reference sequence:\n")
    query = input("Please input the query sequence:\n")
    
    # 输入验证：确保输入的所有字符均为 A/T/G/C
    valid_bases = set("ATGC") # set()函数用于生成集合
    if not set(reference).issubset(valid_bases):
        print("reference输入含有错误碱基. 请只使用大写的 ATGC。")
        return
    if not set(query).issubset(valid_bases):
        print("query输入含有错误碱基. 请只使用大写的 ATGC。")
        return

    # 记录参考序列和查询序列的长度
    m = len(reference)
    n = len(query)
    
    # 一、预处理：构建参考序列所有连续子串的映射字典
    # 定义字典：key为子串，value为一个元组 (起始位置, 是否为互补翻转)
    # 注意：若子串在参考序列正链中出现，将记录正链信息；如果正链没有出现而互补翻转中出现，则记录互补翻转信息
    substr_dict = {}
    
    # 第一部分：正链子串枚举
    # 枚举参考序列中所有连续子串，并记录它们在正链中的起始位置
    for i in range(m): # 默认从 0 开始，一直到 m - 1
        sub = ""
        # 内层循环逐字符构造以 i 为开头的子串
        for j in range(i, m):
            sub += reference[j]
            # 保存正链信息，若该子串已存在，则覆盖无影响（因为正链本来就是优先选择的）
            substr_dict[sub] = (i, False) # 因此，如果是ATGCATGC，索引ATGC时，会找到后面的ATGC（前面的被覆盖）
    
    # 第二部分：互补翻转子串枚举
    # 先计算参考序列的互补翻转序列
    reverse = comp_rev(reference)
    # 枚举互补翻转序列的所有连续子串
    for i in range(m):
        sub = ""
        for j in range(i, m):
            sub += reverse[j]
            l = j - i + 1  # 子串长度
            # 通过公式还原正链上的起始位置（原起始位置 = m - l - i）
            orig_start = m - l - i
            # 如果该子串正链中未出现，则记录互补翻转信息
            if sub not in substr_dict:
                substr_dict[sub] = (orig_start, True)
    
    # 二、动态规划：划分查询序列，寻找最优的切分方案
    # 定义dp数组，其中dp[i]记录从查询序列下标i到末尾的最小切分数量
    # 初始化为 n+1 (一个足够大的值)，并设置dp[n] = 0（空序列需要0次切分）
    dp = [n + 1] * (n + 1)
    dp[n] = 0
    
    # cutting用于记录切分路径：cutting[i] = (下一个切分位置, (参考序列起始位置, 是否互补))
    cutting = [None] * (n + 1) # 列表，长度为 n + 1，每个元素是一个元组
    
    # 倒序遍历查询序列，从末尾向前寻找每个位置的最优切分
    for i in range(n - 1, -1, -1): # 从 n - 1 到 0
        sub = ""
        # 从位置i开始逐步构造子串，并判断是否能在参考序列中找到
        for j in range(i, n):
            sub += query[j]
            if sub in substr_dict: # 找到匹配的子串
                # 如果使用当前子串作为一个切分段能使切分总数更少，则更新dp和cutting信息
                if dp[i] > 1 + dp[j + 1]:
                    dp[i] = 1 + dp[j + 1] # 更新切分总数，使之更少
                    cutting[i] = (j + 1, substr_dict[sub]) # 更新切分路径信息
                    # 这里的j + 1是指下一个切分的位置，substr_dict[sub]则指从当前位置 i 开始，
                    # 找到的最佳后继子串 sub 在参考序列中的信息
    
    # 如果没有找到匹配的切分方案，则输出错误信息并退出
    if cutting[0] is None: # 如果从 0 开始的切分元组为空
        print("Error: 无法匹配query序列中的任何子串，请检查reference和query是否正确。")
        return
    
    # 三、根据cutting数组恢复查询序列的匹配分段信息
    # parts列表用于保存每个匹配分段的信息，格式为((ref_start, ref_end), inversion)
    parts = []
    pos = 0 # 从0开始，往后找
    # 按照cutting记录从查询序列最前面逐段还原匹配结果
    while pos < n:
        nxt, (start, inv) = cutting[pos]
        seg_len = nxt - pos # 当前匹配段在查询序列中的长度
        # 对应在参考序列中的匹配区间为 [start, start + seg_len - 1]
        parts.append(((start, start + seg_len - 1), inv)) # 加入到parts列表中，为后续output做准备
        pos = nxt   # 移到下一个未匹配的位置
    
    # 四、合并连续相同的匹配段：相同的参考起始位置和翻转标志连续出现时合并为一条记录
    result = []
    i = 0
    while i < len(parts):
        # 当前段的参考序列起始位置
        cur_start = parts[i][0][0] # 即start
        # 单位长度: 当前匹配段在参考序列上的长度 (参考结束位置 - 起始位置 + 1)
        unit_len = parts[i][0][1] - parts[i][0][0] + 1 # 即seg_len
        cur_inv = parts[i][1] # 即inv
        count = 1  # 初始化重复次数为1
        j = i + 1 # j表示下一个匹配段的索引
        # 如果后续的连续匹配段与当前段匹配信息相同，则计数累加
        while j < len(parts) and parts[j][0][0] == cur_start and parts[j][1] == cur_inv:
            count += 1
            j += 1
        # 将合并后的结果保存为字典
        result.append({
            'start_pos': cur_start,
            'length': unit_len,
            'repeat_count': count,
            'is_inversion': cur_inv
        })
        i = j  # 更新索引，并跳过已合并的段
    
    # 五、输出
    print("Notice:")
    print("1. 本输出格式中，从前往后的顺序考察query序列中的每一段分别来自于reference的何位置，复制多少次，是否翻转等信息。")
    print("2. 每条记录中，start_pos 表示 query 的该段序列在 reference 里的起始位置（0 为起始索引位置），length 为重复单位长度，repeat_count 为重复次数（至少 1 次，即其本身），is_inversion 表示该段序列是否翻转。") 
    for ans in result:
        print(ans)

if __name__ == "__main__":
    main()