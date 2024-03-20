use assert_cmd::prelude::*; // Add methods on commands
use assert_fs::prelude::*;
use predicates::prelude::*;
use std::process::Command;

#[test]
fn test_missing_arguments() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("sumdiff")?;
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Usage: sumdiff <left> <right>"));
    Ok(())
}
#[test]
fn test_missing_files() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("sumdiff")?;
    cmd.arg("test/a").arg("test/b");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Unable to open test/a"));
    Ok(())
}

#[test]
fn test_null_diff() -> Result<(), Box<dyn std::error::Error>> {
    let a_file = assert_fs::NamedTempFile::new("a.txt")?;
    a_file.write_str("Hello world\n")?;
    let b_file = assert_fs::NamedTempFile::new("b.txt")?;
    b_file.write_str("Hello world\n")?;
    let mut cmd = Command::cargo_bin("sumdiff")?;
    cmd.arg(a_file.path()).arg(b_file.path());
    cmd.assert().success().stdout(predicate::str::is_empty());
    cmd.assert().success().stderr(predicate::str::is_empty());
    Ok(())
}


#[test]
fn test_single_diff() -> Result<(), Box<dyn std::error::Error>> {
    let a_file = assert_fs::NamedTempFile::new("a.txt")?;
    a_file.write_str("Goodbye\ncruel\nworld\n")?;
    let b_file = assert_fs::NamedTempFile::new("b.txt")?;
    b_file.write_str("Goodbye\nsweet\nworld\n")?;
    let mut cmd = Command::cargo_bin("sumdiff")?;
    cmd.arg(a_file.path()).arg(b_file.path());
    cmd.assert().success().stdout(predicate::str::ends_with("\
@@ -1,3 +1,3 @@
 Goodbye
-cruel
+sweet
 world
"));
    cmd.assert().success().stderr(predicate::str::is_empty());
    Ok(())
}

#[test]
fn test_missing_newline_new() -> Result<(), Box<dyn std::error::Error>> {
    let a_file = assert_fs::NamedTempFile::new("a.txt")?;
    a_file.write_str("Goodbye\ncruel\nworld\n")?;
    let b_file = assert_fs::NamedTempFile::new("b.txt")?;
    b_file.write_str("Goodbye\nsweet\nworld")?;
    let mut cmd = Command::cargo_bin("sumdiff")?;
    cmd.arg(a_file.path()).arg(b_file.path());
    cmd.assert().success().stdout(predicate::str::ends_with("\
@@ -1,3 +1,3 @@
 Goodbye
-cruel
-world
+sweet
+world
\\ No newline at end of file
"));
    cmd.assert().success().stderr(predicate::str::is_empty());
    Ok(())
}


#[test]
fn test_missing_newline_old() -> Result<(), Box<dyn std::error::Error>> {
    let a_file = assert_fs::NamedTempFile::new("a.txt")?;
    a_file.write_str("Goodbye\ncruel\nworld")?;
    let b_file = assert_fs::NamedTempFile::new("b.txt")?;
    b_file.write_str("Goodbye\nsweet\nworld\n")?;
    let mut cmd = Command::cargo_bin("sumdiff")?;
    cmd.arg(a_file.path()).arg(b_file.path());
    cmd.assert().success().stdout(predicate::str::ends_with("\
@@ -1,3 +1,3 @@
 Goodbye
-cruel
-world
\\ No newline at end of file
+sweet
+world
"));
    cmd.assert().success().stderr(predicate::str::is_empty());
    Ok(())
}