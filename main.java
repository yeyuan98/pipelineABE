
/**
 * Run rABE analysis of a single sample from alignment BAM to variant VCF
 *   toolchain idea originally from https://github.com/linyz/remora
 *
 * @author Ye Yuan
 * @version 0
 */

import java.lang.ProcessBuilder;
import java.nio.file.Path;
import java.nio.file.Files;
import java.io.File;
import java.io.IOException;
import java.nio.file.FileAlreadyExistsException;
import java.util.Arrays;
import java.util.List;

import java.time.Duration;
import java.time.Instant;

public class main
{
    private static void checkArguments(String[] args){
        if ( args.length != 3 ){
            System.out.println("Usage: pipelineABE <type = (dna, rna)> <Alignment BAM> <Output Directory>");
            System.exit(1);
        }
        
        Path bamPath = Path.of(args[1]);
        Path outDir = Path.of(args[2]);
        
        if ( !Files.isReadable(bamPath) ){
            System.out.println("Sample BAM is not readable.");
            System.exit(2);
        }
        
        try {
            Files.createDirectory(outDir);
        }
        catch( FileAlreadyExistsException e ){
            System.out.println("Output directory already exists.");
            //System.exit(2);
        }
        catch( IOException e){
            System.out.println("IO error while creating output directory.");
            System.exit(2);
        }
    }
    
    public static void main(final String[] args) throws InterruptedException{
        /*
         * args[] = (BAM path, output_dir)
         */
        // Check arguments and create output dir
        Instant timeStart = Instant.now();
        Instant timeNow;
        checkArguments(args);
        String type = args[0];
        Path bamPath = Path.of(args[1]);
        Path outDir = Path.of(args[2]);
        
        // samtools index
        timeNow = Instant.now();
        System.out.println(
        String.format("[%s] Indexing input BAM...", 
        Duration.between(timeStart, timeNow)));
        
        Path samtoolsBinPath = Path.of("/Users/yeyuan/mambaforge/envs/binfo/bin/samtools");
        ProcessBuilder pb = 
            new ProcessBuilder(samtoolsBinPath.toString(), "index", "-@ 6", bamPath.toString());
        try{
            Process p = pb.start();
            p.waitFor();
        } catch ( IOException e ) {
            System.out.println("IO error while running samtools index");
        }
        
        // picard MarkDuplicate
        timeNow = Instant.now();
        System.out.println(
        String.format("[%s] Marking duplicate...", 
        Duration.between(timeStart, timeNow)));
        
        Path JavaPath = Path.of("/usr/bin/java");
        Path picardJarPath = Path.of("/Users/yeyuan/mambaforge/envs/binfo/share/picard-3.0.0-1/picard.jar");
        Path picardLogPath = outDir.resolve("MarkDuplicates.log");
        Path picardOutPath = outDir.resolve("process.MarkDuplicate.bam");
        Path picardMetricPath = outDir.resolve("MarkDuplicates.metrics.txt");
        pb = 
            new ProcessBuilder(JavaPath.toString(), 
            "-jar", picardJarPath.toString(), "MarkDuplicates", "-I", bamPath.toString(), 
            "-M", picardMetricPath.toString(), "-O", picardOutPath.toString());
        pb.redirectError(ProcessBuilder.Redirect.to(
                picardLogPath.toFile()
            ));
        try{
            Process p = pb.start();
            p.waitFor();
            
            if (p.exitValue() != 0){
                System.out.println("Error occured while marking duplicates. Check log.");
                System.exit(1);
            }
        } catch ( IOException e ) {
            System.out.println("IO error while marking duplicates");
            e.printStackTrace();
        }
        
        // GATK SplitNCigarReads - only for RNA sequencing
        Path gatkOutPath = picardOutPath;
        Path genomeRefPath = Path.of("/Volumes/SeagateExpansionYY/Ye/sequencing/RNAseq/enc_rABE/index/dm6.fa");
        if (type.equals("rna")){
            
            timeNow = Instant.now();
            System.out.println(
            String.format("[%s] Splitting reads spanning intron...", 
            Duration.between(timeStart, timeNow)));
            
            Path gatkJarPath = Path.of("/Users/yeyuan/mambaforge/envs/binfo/share/gatk4-4.4.0.0-0/gatk-package-4.4.0.0-local.jar");
            Path gatkLogPath = outDir.resolve("SplitNCigarReads.log");
            gatkOutPath = outDir.resolve("process.SplitNCigarReads.bam");
            pb = 
                new ProcessBuilder(JavaPath.toString(),
                "-jar", gatkJarPath.toString(), "SplitNCigarReads",
                "-I", picardOutPath.toString(), "-O", gatkOutPath.toString(),
                "-R", genomeRefPath.toString());
            pb.redirectError(ProcessBuilder.Redirect.to(
                    gatkLogPath.toFile()
                ));
            try{
                Process p = pb.start();
                p.waitFor();
                
                if (p.exitValue() != 0){
                    System.out.println("Error occured while splitting intron reads. Check log.");
                    System.exit(1);
                }
            } catch ( IOException e ) {
                System.out.println("IO error while splitting intron reads");
                e.printStackTrace();
            }
        }
    
        // bcftools variant calling
        timeNow = Instant.now();
            System.out.println(
            String.format("[%s] Calling variants against reference genome...", 
            Duration.between(timeStart, timeNow)));
        
        Path bcftoolsBinPath = Path.of("/Users/yeyuan/mambaforge/envs/binfo/bin/bcftools");
        Path bcftoolsOutPath = outDir.resolve("mpileup.vcf");
        Path bcftoolsLogPath = outDir.resolve("bcftools.log");
        String bcftoolsAnnotations = "DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR";
        String bcftoolsFilter = "INFO/AD[1-]>2 & MAX(FORMAT/DP)>20";
        
        ProcessBuilder[] pbs = {
            new ProcessBuilder(
            bcftoolsBinPath.toString(), "mpileup", 
            "-f", genomeRefPath.toString(),
            "-d", "10000000", "-I", "-a", bcftoolsAnnotations,
            gatkOutPath.toString(), "-O", "v"),
            new ProcessBuilder(
            bcftoolsBinPath.toString(), "filter",
            "-i", bcftoolsFilter, "-O", "v", "-")
        };
        
        pbs[1].redirectOutput(bcftoolsOutPath.toFile());
        pbs[0].redirectError(bcftoolsLogPath.toFile()); // filter yields no log output in error stream
        try{
            
            List<Process> pss = ProcessBuilder.startPipeline(Arrays.asList(pbs));
            pss.get(1).waitFor();
            
            if (pss.get(1).exitValue() != 0){
                System.out.println("Error occured while calling variants.");
                System.exit(1);
            }
        } catch ( IOException e ){
            System.out.println("IO error while running bcftools");
                e.printStackTrace();
        }
        
        // done
        timeNow = Instant.now();
            System.out.println(
            String.format("[%s] Done.", 
            Duration.between(timeStart, timeNow)));
    }
    
}
