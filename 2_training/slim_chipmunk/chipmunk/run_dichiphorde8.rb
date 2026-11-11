#!/usr/bin/ruby
$: << "#{File.dirname(File.expand_path(__FILE__))}/ytilib/"
require "ytilib.rb"

ns = ['A','C','G','T']
$dins = ns.collect { |n| ns.collect { |n1| "#{n}#{n1}" } }.flatten  

class PM
  def rcifn!
    rcdone = false
    g,c = @matrix['G'].sum, @matrix['C'].sum
    a,t = @matrix['A'].sum, @matrix['T'].sum
    if g+c > a+t
      (dummy = self.revcomp!, rcdone = true) if g <= c
    else
      (dummy = self.revcomp!, rcdone = true) if t <= a
    end
    return rcdone    
  end
end

def chipstep(motif_name_t, params)

  cmdline = "java -Xms2g -Xmx10g -cp \"#{File.dirname(File.expand_path(__FILE__))}/chipmunk.jar\" ru.autosome.di.ChIPHorde #{params}"
  report "running dichipmunk-engine: #{cmdline}"
  
  result_full = `#{cmdline}`.strip
  File.open("#{motif_name_t}_dichipmunk.log", "w") { |f|
    f << result_full
  }
  
  $diagnosis = result_full.split("INFO|")[1].split("\n")[0]
  
  result_full = result_full.split("OUTC|ru.autosome.di.ChIPHorde")[1].strip
  result_full = result_full.split("MOTF")
  result_full.shift

  result_full.each_with_index { |result, ri|
  
    motif_name = "#{motif_name_t}_#{ri}"
    
    result = result.strip.split("\n")
    result_h = {}
    
    words = {}
    result.each { |line|
      key,value = line.split("|")
      next if key != "WORD"
      length = value.split[2].size
      words[length] = [] unless words[length]
      words[length] << value
    }
    
    result.each { |v|
      key, value = v.split("|")
      if result_h[key]
        if result_h[key].is_a?(Array)
          result_h[key] << value
        else
          result_h[key] = [result_h[key], value]
        end
      else
        result_h[key] = value
      end
    }
    
    #$diagnosis = result_h["DIAG"]
    words_count = result_h["N"].to_f
    
    File.open("#{motif_name}.dpcm", "w") { |f|
      $dins.each { |din|
        f.puts("#{result_h[din]}")
      }
    }
    File.open("#{motif_name}.dpwm", "w") { |f|
      $dins.each { |din|
        pw = "PW#{din}"
        f.puts("#{result_h[pw]}")
      }
    }
    
    File.open("#{motif_name}_words.txt", "w") { |f|
      words[words.keys.last].each { |w| 
        f.puts(w)
      }
    } unless words.empty?
    
    File.open("#{motif_name}.dat","w") { |f|
      f.puts("TOTAL = #{result_h['TOTL']}")
      f.puts("SEQUENCES = #{result_h['SEQS']}")
      f.puts("WORDS = #{result_h['WRDS']}")
      f.puts("KDIDIC = #{result_h['KDDC']}")
      f.puts("THRESHOLD = #{result_h['THRE']}")
      f.puts("PVALUE = #{result_h['PVAL']}")
      diag = [result_h['DIAG']].flatten.first
      f.puts("DIAGNOSIS = #{diag}")
    }
    
    system("ruby #{File.dirname(File.expand_path(__FILE__))}/dpmflogo3.rb #{motif_name}.dpcm #{motif_name}.png")
  
  }
end

report "run_dichiphorde7.rb started, usage <motif_name> <ChIPHorde-engine-params>"
start __FILE__

exit(2) if ARGV.size < 2

motif_name = ARGV.shift

$rna = ARGV.inject(false) { |sum,v| sum |= v[0,2].downcase == "r:" }

params = ARGV.inject("") { |pa, pr| pa += "#{pr} " }

chipstep(motif_name, params)

report $diagnosis
