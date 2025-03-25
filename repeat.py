def comp_rev(seq):
    """计算DNA序列的互补翻转"""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(comp[c] for c in reversed(seq))

def main():
    ref_seq = input("Please input the reference sequence:\n")
    qry_seq = input("Please input the query sequence:\n")
    
    m = len(ref_seq)
    n = len(qry_seq)
    # 构建字典：子串 -> (起始位置, 是否翻转)
    substr_dict = {}
    for i in range(m):
        sub = ""
        for j in range(i, m):
            sub += ref_seq[j]
            substr_dict[sub] = (i, False)
    
    ref_rc = comp_rev(ref_seq)
    for i in range(m):
        sub = ""
        for j in range(i, m):
            sub += ref_rc[j]
            l = j - i + 1
            orig_start = m - l - i
            if sub not in substr_dict:
                substr_dict[sub] = (orig_start, True)
    
    # 动态规划划分 qry_seq
    dp = [n+1]*(n+1)
    dp[n] = 0
    trace = [None]*(n+1)  # trace[i] = (next_index, (start_pos, is_inversion))
    
    for i in range(n-1, -1, -1):
        sub = ""
        for j in range(i, n):
            sub += qry_seq[j]
            if sub in substr_dict:
                if dp[i] > 1 + dp[j+1]:
                    dp[i] = 1 + dp[j+1]
                    trace[i] = (j+1, substr_dict[sub])
    
    # 根据 trace 得到分段信息，每段为 (reference 起始, reference 结束, is_inversion)
    segments = []
    pos = 0
    while pos < n:
        nxt, (start, inv) = trace[pos]
        seg_len = nxt - pos
        segments.append(((start, start + seg_len - 1), inv))
        pos = nxt

    # 合并连续相同(起始位置及翻转标志相同)的段，统计重复次数
    outputs = []
    i = 0
    while i < len(segments):
        cur_start = segments[i][0][0]
        unit_len = segments[i][0][1] - segments[i][0][0] + 1
        cur_inv = segments[i][1]
        count = 1
        j = i + 1
        while j < len(segments) and segments[j][0][0] == cur_start and segments[j][1] == cur_inv:
            count += 1
            j += 1
        outputs.append({
            'start_pos': cur_start,
            'length': unit_len,
            'repeat_count': count,
            'is_inversion': cur_inv
        })
        i = j

    print("Notice:")
    print("1. 本输出格式中，从前往后的顺序考察query序列中的每一段分别来自于reference的何位置，复制多少次，是否翻转等信息。")
    print("2. 每条记录中，start_pos 表示 query 的该段序列在 reference 里的起始位置（0 为起始索引位置），length 为重复单位长度，repeat_count 为重复次数（至少 1 次，即其本身），is_inversion 表示该段序列是否翻转。")
    
    for out in outputs:
        print(out)

if __name__ == "__main__":
    main()